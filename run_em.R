##############################################
## run_em function (EM algoritm)            ##
## Density Estimation                       ##
## author: karolina.kulec@ucdconnect.ie     ##
## date: June 2019                          ##
##############################################

source("auxillary_functions.R")

run_em <- function(data,
                   groups,
                   max_iterations = 200,
                   convergence_limit = 1e-05,
                   random_initialisation_runs = 50) {
  # Runs the EM algorithm for a finite mixture of unimodal beta distributions.
  #
  # Args:
  #   data: data vector believed to originate from a mixture of unimodal betas.
  #   groups: number of groups/components.
  #   max_iterations: maximum number of algorithm iterations.
  #   convergence_limit: lack of progress measure indicating convergence of the algorithm.
  #   random_initialisation_runs: the number of random runs to be passed to function initialise_me.
  # Returns:
  #   output: output list of class 'BetaEM'
  
  N <- length(data)
  log_likelihood <- NULL
  initial_values <-
    initialise_me(data = data,
                  random_runs = random_initialisation_runs,
                  groups = groups)
  
  if (groups == 1) {
    theta <- initial_values
    log_likelihood[1] <-
      sum(log(dbeta.rep(data, theta[1], theta[2])))
    optim_theta <-
      c(inv_location_constraint(theta[1]), log(theta[2]))
    maximise <-
      optim(optim_theta,
            fn = max_me,
            data = data,
            groups = groups)
    if (maximise$convergence == 1) {
      stop("optim: iteration limit maxit had been reached without convergence")
    }
    theta <-
      c(location_constraint(maximise$par[1]), exp(maximise$par[2]))
    log_likelihood[2] <-
      sum(log(dbeta.rep(data, theta[1], theta[2])))
    
    iterations <- 2  # 1 + initial values
    
  } else {
    theta <- array(0, c(groups, 2, max_iterations + 1))
    mixing_proportions  <- array(0, c(groups, max_iterations + 1))
    by_group <- NULL
    z <- NULL
    
    theta[, , 1] <- initial_values[[1]]
    mixing_proportions[, 1] <- initial_values[[2]]
    
    log_likelihood[1] <- sum(log(
      sapply(1:groups, function(j)
        dbeta.rep(data, location = theta[j, 1, 1], scale = theta[j, 2, 1])) %*% mixing_proportions[, 1]
    ))
    
    by_group  <-
      sapply(1:groups, function(j)
        mixing_proportions[j, 1] %*% dbeta.rep(data, location = theta[j, 1, 1], scale = theta[j, 2, 1]))
    prior_mix <- apply(by_group, 1, sum)
    z <- by_group / rep(prior_mix, groups)
    
    iterations <- 1
    
    for (r in 1:max_iterations) {
      # E step
      last_z <- z
      prior_mix <- apply(by_group, 1, sum)
      z <- by_group / rep(prior_mix, groups)
      
      # M step
      # convert the parameter matrix to vector for optim compatability
      optim_theta <-
        cbind(inv_location_constraint(theta[, 1, r]), log(theta[, 2, r]))
      optim_vec <- numeric()
      for (i in 1:groups) {
        optim_vec <- c(optim_vec, optim_theta[i,])
      }
      
      maximise <-
        optim(
          par = optim_vec,
          fn = max_me,
          method = "BFGS",
          data = data,
          z = z,
          groups = groups
        )
      new_par <- maximise$par
      if (maximise$convergence == 1) {
        stop("optim: iteration limit maxit had been reached without convergence")
      }
      
      # convert the parameter vector to parameter matrix
      new_params <- matrix(new_par, nrow = groups, byrow = T)
      new_theta <-
        cbind(location_constraint(new_params[, 1]), exp(new_params[, 2]))
      new_mixing_proportions <- apply(z, 2, sum) / N
      
      # Updating step
      theta[, , (r + 1)] <- new_theta
      mixing_proportions[, (r + 1)] <- new_mixing_proportions
      
      log_likelihood[r + 1] <- sum(log(
        sapply(1:groups, function(j)
          dbeta.rep(data, location = theta[j, 1, r + 1], scale = theta[j, 2, r + 1])) %*% mixing_proportions[, r + 1]
      ))
      
      by_group <- NULL
      by_group <- sapply(1:groups, function(j)
        mixing_proportions[j, r + 1] %*% dbeta.rep(data, location = theta[j, 1, r + 1], scale = theta[j, 2, r + 1]))
      
      
      if (abs(log_likelihood[r + 1] - log_likelihood[r]) < convergence_limit) {
        break
      }
      iterations <- iterations + 1
      if (iterations == max_iterations + 1) {
        warning(
          paste(
            "Max iteration limit",
            max_iterations,
            "has been reached without convergence"
          )
        )
      }
    }
  }
  
  AIC <- -2 * log_likelihood[iterations] + 2 * (3 * groups - 1)
  BIC <- -2 * log_likelihood[iterations] + (3 * groups - 1) * log(N)
  
  if (groups == 1) {
    output <- list(
      "Theta" = theta,
      "AIC" = AIC,
      "BIC" = BIC,
      "LogLikelihood" = log_likelihood,
      "Data" = data,
      "Groups" = groups
    )
  } else {
    output <- list(
      "Iterations" = iterations,
      "Theta" = theta[, , iterations],
      "MixingProportions" = mixing_proportions[, iterations],
      "AIC" = AIC,
      "BIC" = BIC,
      "LogLikelihood" = log_likelihood,
      "ClassProbabilities" = z,
      "Data" = data,
      "Groups" = groups
    )
  }
  
  cat("Object of class 'BetaEM' \n")
  cat("======================== \n")
  cat("Available components \n")
  print(names(output))
  
  class(output) <- "BetaEM"
  return(output)
}


