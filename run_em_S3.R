##############################################
## S3 methods for run_em                    ##
##                                          ##
##  1. summary                              ##  
##  2. plot                                 ##
##  3. label_switch                         ##
##  4. beta_clus                            ##
##                                          ##
## author: karolina.kulec@ucdconnect.ie     ##
## date: 10 June 2019                       ##
##############################################

source("run_em.R")

summary.BetaEM <- function(x) {
  # Summary of run_em output.
  #
  # Args:
  #   x: object of class BetaEM.
  # Returns:
  # Number of iterations until convergence.
  # Final mixture parameters (location, scale and mixin proportions).
  stopifnot(inherits(x, "BetaEM"))
  
  location <- x$Theta[, 1]
  scale    <- x$Theta[, 2]
  
  cat('EM iterations until convergence \n')
  print(x$Iterations)
  cat(sep = '\n')
  cat('-- Final mixture parameters  -- \n', sep = "\n")
  cat("Location \n")
  print(location)
  cat(sep = "\n")
  cat("Scale \n")
  print(scale)
  cat(sep = "\n")
  cat("Mixing Proportions \n")
  print(x$MixingProportions)
}


plot.BetaEM <- function(x, groups = 3) {
  # Diagnostic plots for run_em.
  #
  # Args:
  #   x: object of class BetaEM.
  #   data: simulated data or other data on which the algorithm was used
  #   groups: number of components/groups in the mixture.
  # Returns:
  #   Histogram of data with estimated density overlay.
  #   Plot of the log likelihood vs. number of iterations.
  stopifnot(inherits(x, "BetaEM"))
  
  location           <- x$Theta[, 1]
  scale              <- x$Theta[, 2]
  
  evaluation_points  <- seq(0, 1, length.out = 1000)
  estimated_density  <-
    sapply(1:groups, function(j)
      dbeta.rep(evaluation_points, location[j], scale[j])) %*% x$MixingProportions
  
  par(mfrow = c(1, 2))
  
  # plot 1
  hist(
    x$Data,
    breaks = 25,
    freq = F,
    xlim = c(0, 1),
    xlab = "",
    main = "Data histogram with estimated density overlay",
    cex.main = 0.9
  )
  lines(evaluation_points, estimated_density, col = "red")
  
  # plot 2
  plot(
    x$LogLikelihood,
    type = "b",
    xlab = "Iterations",
    ylab = "Log Likelihood",
    main = "Log likelihood vs. number of iterations",
    cex.main = 0.9
  )
}

label_switch <-
  function (x,
            SimulationComponents,
            groups = 3,
            N = 1000) {
    UseMethod("label_switch", x)
  }

label_switch.BetaEM <-
  function(x,
           SimulationComponents,
           groups = 3,
           N = 1000) {
    # Solves the label switching problem. See http://www-personal.k-state.edu/~wxyao/material/submitted/labelswitchingfrequency.pdf.
    #### Flagged for problems ####
    # Args:
    #   x: object of class BetaEM.
    #   SimulationComponents: simulated data components.
    #   groups: number of components/groups in the mixture.
    #   N: length of data.
    # Returns:
    #   true_labels: correct label vector.
    stopifnot(inherits(x, "BetaEM"))
    
    permutation_matrix <- permute_me(groups)
    permutations <- nrow(permutation_matrix)
    
    z <- matrix(0, nrow = N, ncol = groups)
    z <-
      sapply(1:groups, function(i)
        z[, i] <- ifelse(SimulationComponents == i, 1, 0))
    
    class_probabilities <- x$ClassProbabilities #(N x groups matrix)
    permute_class_probabilities  <-
      array(0, c(N, groups, permutations))
    for (i in 1:permutations) {
      permute_class_probabilities[, , i] <-
        class_probabilities[, permutation_matrix[i, ]]
    }
    
    maximise <-
      sapply(1:permutations, function(i)
        sum(sapply(1:groups, function(j)
          z[, j] * log(permute_class_probabilities[, j, i]))))
    
    true_labels <-
      as.vector(permutation_matrix[which.max(maximise),])
    
    # switching up the problem labels...
    if (all(true_labels == c(2, 3, 1))) {
      true_labels <- c(3, 1, 2)
    } else if (all(true_labels == c(3, 1, 2))) {
      true_labels <- c(2, 3, 1)
    } else {
      true_labels
    }
    
    return(true_labels)
  }


beta_clus <-
  function (x,
            SimulationComponents,
            PlotUncertainties = TRUE) {
    UseMethod("beta_clus", x)
  }

beta_clus.BetaEM <- function(x, SimulationComponents) {
  # Gives the clustering solution based on run_em.
  #
  # Args:
  #   x: object of class BetaEM.
  #   SimulationComponents: simulated data components.
  # Returns:
  #   Cross tabulation of simulated data components and solution componenets (class memberships).
  #   Adjusted Rand Index (ARI).
  #   Class memberships.
  #   Class membership uncertainties.
  true_labels <- label_switch(x, SimulationComponents)
  class_probabilities <- x$ClassProbabilities
  class_membership_uncertainties <-
    apply(class_probabilities, 1, function(x)
      1 - max(x))
  class_memberships <- apply(class_probabilities, 1, which.max)
  relabelled_class_memberships <-
    ifelse(
      class_memberships == 1,
      true_labels[1],
      ifelse(class_memberships == 2, true_labels[2], true_labels[3])
    )
  cross_tab <-
    table(" " <-
            SimulationComponents,
          " " <- relabelled_class_memberships)
  ARI <- adjusted_rand(cross_tab)
  
  list(
    "CrossTab" = cross_tab,
    "AdjustedRandIndex" = ARI,
    "ClassMemberships" = relabelled_class_memberships,
    "ClassMembershipUncertainties" = class_membership_uncertainties
  )
  
}
