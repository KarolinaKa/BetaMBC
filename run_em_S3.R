##############################################
## S3 methods for run_em                    ##
##                                          ##
##  1. summary                              ##
##  2. plot                                 ##
##  3. label_switch                         ##
##  4. beta_clus                            ##
##  5. simulation_beta_clus                 ##
##                                          ##
## author: karolina.kulec@ucdconnect.ie     ##
## date: 10 June 2019                       ##
##############################################

source('run_em.R')


summary.BetaEM <- function(x) {
  # Summary of run_em output.
  #
  # Args:
  #   x: object of class BetaEM.
  # Returns:
  # Summary of BetaEM object.
  stopifnot(inherits(x, 'BetaEM'))
  if (x$Groups == 1) {
    cat('-- Final density parameters  -- \n', sep = '\n')
    cat('Location \n')
    cat(x$Theta[1], '\n')
    cat(sep = '\n')
    cat('Scale \n')
    cat(x$Theta[2], '\n')
    
  } else {
    cat(paste('EM iterations until convergence:', x$Iterations),
        sep = '\n')
    
    cat(sep = '\n')
    cat('-- Final mixture parameters  -- \n', sep = '\n')
    cat('Location \n')
    cat(x$Theta[, 1], '\n')
    cat(sep = '\n')
    cat('Scale \n')
    cat(x$Theta[, 2], '\n')
    cat(sep = '\n')
    cat('Mixing Proportions \n')
    cat(x$MixingProportions, '\n')
  }
  cat(sep = '\n')
  cat(paste('Number of component densities:', x$Groups), sep = '\n')
  cat(paste('Number of estimated parameters:', 3 * x$Groups - 1), sep = '\n')
  cat(paste('Maximum log likelihood:', max(x$LogLikelihood)), sep = '\n')
  cat(sep = '\n')
  
  cat(paste('BIC:', x$BIC), sep = '\n')
  cat(paste('AIC:', x$AIC), sep = '\n')
}


plot.BetaEM <- function(x, groups = x$Groups) {
  # Diagnostic plots for run_em.
  #
  # Args:
  #   x: object of class BetaEM.
  #   data: simulated data or other data on which the algorithm was used.
  #   groups: number of components/groups in the mixture.
  # Returns:
  #   Histogram of data with estimated density overlay.
  #   Plot of the log likelihood vs. number of iterations.
  stopifnot(inherits(x, 'BetaEM'))
  
  evaluation_points  <- seq(0, 1, length.out = 1000)
  
  if (groups == 1) {
    estimated_density <-
      dbeta.rep(evaluation_points, x$Theta[1], x$Theta[2])
    
  } else {
    location           <- x$Theta[, 1]
    scale              <- x$Theta[, 2]
    
    estimated_density  <-
      sapply(1:groups, function(j)
        dbeta.rep(evaluation_points, location[j], scale[j])) %*% x$MixingProportions
  }
  
  par(mfrow = c(1, 2))
  
  # plot 1
  hist(
    x$Data,
    breaks = 25,
    freq = F,
    xlim = c(0, 1),
    xlab = '',
    main = 'Data histogram with estimated density overlay',
    cex.main = 0.9
  )
  lines(evaluation_points,
        estimated_density,
        type = 'l',
        col = 'red', 
        lwd = 2)
  
  # plot 2
  plot(
    x$LogLikelihood,
    type = 'b',
    xlab = 'Iterations',
    ylab = 'Log Likelihood',
    main = 'Log likelihood vs. number of iterations',
    cex.main = 0.9
  )
  
  par(mfrow = c(1, 1))

}

simulation_plot <-
  function(x, 
           simulation) {
    UseMethod("simulation_plot", x)
  }

simulation_plot <-
  function(x,
           simulation) {
    stopifnot(inherits(x, 'BetaEM'))
    stopifnot(inherits(simulation, "Simulation"))
    # Compares the true and estimated densities.
    #
    # Args:
    #   x: object of class BetaEM.
    #   simulation: object of class Simulation. 
    # Returns:
    #   Histogram of simulated data with estimated and true density overlay.
    evaluation_points  <- seq(0, 1, length.out = 1000)
    
    location <- x$Theta[, 1]
    scale <- x$Theta[, 2]
    
    estimated_density  <-
      sapply(1:x$Groups, function(j)
        dbeta.rep(evaluation_points, location[j], scale[j])) %*% x$MixingProportions
    component1 <-
      dbeta.rep(evaluation_points, location[1], scale[1])
    true_density <-
      sapply(1:x$Groups, function(j)
        dbeta.rep(evaluation_points, simulation$Location[j], simulation$Scale[j])) %*% simulation$MixingProportions
    hist(
      simulation$SimulatedData,
      breaks = 25,
      freq = F,
      xlim = c(0, 1),
      xlab = '',
      main = 'Data histogram with true and estimated density overlay',
      cex.main = 0.9
    )
    lines(
      evaluation_points,
      estimated_density,
      type = 'l',
      col = 'red',
      lwd = 2
    )
    lines(
      evaluation_points,
      true_density,
      type = "l",
      lwd = 2,
      lty = 2
    )
    legend(
      "topright",
      col = c("red", "black"),
      lty = 1:2,
      lwd = 2,
      legend = c("Estimate", "True"),
      cex = 0.8,
      bty = "n"
    )
  }

label_switch <-
  function (x,
            simulation)
  {
    UseMethod('label_switch', x)
  }

label_switch.BetaEM <-
  function(x,
           simulation)
  {
    # Solves the label switching problem. See http://www-personal.k-state.edu/~wxyao/material/submitted/labelswitchingfrequency.pdf.
    #### Flagged for problems ####
    # Args:
    #   x: object of class "BetaEM".
    #   simulation: object of class "Simulation".
    # Returns:
    #   true_labels: correct label vector.
    stopifnot(inherits(x, 'BetaEM'))
    stopifnot(inherits(simulation, "Simulation"))
    
    groups <- x$Groups
    
    if (groups == 1) {
      stop('label_switch not applicable.')
    }
    
    N <- length(x$Data)
    
    permutation_matrix <- permute_me(groups)
    permutations <- nrow(permutation_matrix)
    
    z <- matrix(0, nrow = N, ncol = groups)
    z <-
      sapply(1:groups, function(i)
        z[, i] <- ifelse(simulation$Components == i, 1, 0))
    
    class_probabilities <- x$ClassProbabilities #(N x groups matrix)
    permute_class_probabilities  <-
      array(0, c(N, groups, permutations))
    for (i in 1:permutations) {
      permute_class_probabilities[, , i] <-
        class_probabilities[, permutation_matrix[i,]]
    }
    
    maximise <-
      sapply(1:permutations, function(i)
        sum(sapply(1:groups, function(j)
          z[, j] * log(permute_class_probabilities[, j, i]))))
    
    true_labels <-
      as.vector(permutation_matrix[which.max(maximise),])
    
    # switching up the problem labels... (only found when groups = 3)
    if (groups == 3) {
      if (all(true_labels == c(2, 3, 1))) {
        true_labels <- c(3, 1, 2)
      } else if (all(true_labels == c(3, 1, 2))) {
        true_labels <- c(2, 3, 1)
      } else {
        true_labels
      }
    }
    
    return(true_labels)
  }


beta_clus <-
  function (x, plot = TRUE) {
    UseMethod('beta_clus', x)
  }

beta_clus.BetaEM <- function(x, plot = TRUE) {
  # Gives the clustering solution based on run_em.
  #
  # Args:
  #   x: object of class BetaEM.
  # Returns:
  #   Class memberships.
  #   Class membership uncertainties.
  #   Plot of class membership uncertainties
  stopifnot(inherits(x, 'BetaEM'))
  if (x$Groups == 1) {
    stop('beta_clus not applicable.')
  }
  class_probabilities <- x$ClassProbabilities
  class_membership_uncertainties <-
    cat("Object of class 'BetaEM' \n")
  cat("======================== \n")
  cat("Available components \n")
  print(names(output))
  
  class(output) <- "BetaEM"
  return(output)
  class_memberships <- apply(class_probabilities, 1, which.max)
  
  if (plot == TRUE) {
    plot(
      class_membership_uncertainties ~ x$Data,
      type = 'h',
      col = adjustcolor("black", alpha.f = 0.5),
      xlab = '',
      xlim = c(0, 1), 
      ylim = c(0, 1),
      ylab = 'Uncertainty',
      main = 'Cluster membership uncertainties (black lines) with estimated density overlay',
      cex.main = 0.9, 
      sub = 'Note: Estimated density is scaled to have unit maximum.',
      cex.sub = 0.8
    )
    
    evaluation_points  <- seq(0, 1, length.out = 1000)
    location           <- x$Theta[, 1]
    scale              <- x$Theta[, 2]
    estimated_density  <-
      sapply(1:x$Groups, function(j)
        dbeta.rep(evaluation_points, location[j], scale[j])) %*% x$MixingProportions
    lines(evaluation_points, estimated_density / max(estimated_density), col = 'red')
  }
  
  list('ClassMemberships' = class_memberships,
       'ClassMembershipUncertainties' = class_membership_uncertainties)
}

simulation_beta_clus <-
  function (x,
            simulation) {
    UseMethod('simulation_beta_clus', x)
  }

simulation_beta_clus.BetaEM <- function(x, simulation) {
  # Gives the clustering solution based on run_em and compares it to the original simulation components.
  #
  # Args:
  #   x: object of class "BetaEM".
  #   simulation: object of class "Simulation".
  # Returns:
  #   Cross tabulation of simulated data components and solution componenets (class memberships).
  #   Adjusted Rand Index (ARI).
  #   Class memberships.
  stopifnot(inherits(x, 'BetaEM'))
  stopifnot(inherits(simulation, "Simulation"))
  
  groups <- x$Groups
  
  if (groups == 1) {
    stop('beta_clus not applicable.')
  }
  
  true_labels <- label_switch(x, simulation)
  class_memberships <- beta_clus(x)$ClassMemberships

  relabelled_class_memberships <-
    ifelse(
      class_memberships == 1,
      true_labels[1],
      ifelse(class_memberships == 2, true_labels[2], true_labels[3])
    )
  cross_tab <-
    table(' ' <-
            simulation$Components,
          ' ' <- relabelled_class_memberships)
  ARI <- adjusted_rand(cross_tab)
  
  list(
    'CrossTab' = cross_tab,
    'AdjustedRandIndex' = ARI,
    'ClassMemberships' = relabelled_class_memberships
  )
  
}
