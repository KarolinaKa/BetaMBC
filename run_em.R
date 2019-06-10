##############################################
## run_em function (EM algoritm)            ##
## author: karolina.kulec@ucdconnect.ie     ##
## date: June 2019                          ##
##############################################

source("auxillary_functions.R")

# Simulate some data 
simulation <- simulate_beta(
  N = 1000,
  mixing_proportions = c(0.5, 0.2, 0.3),
  location = c(0.1, 0.5, 0.9),
  scale = c(0.03, 0.03, 0.04)
)

run_em <- function(data,
                   groups = 3,
                   max_iterations = 200,
                   convergence_limit = 1e-05,
                   random_initialisation_runs = 50) {
  # Runs the EM algorithm for a finite mixture of unimodal beta distributions.
  #
  # Args: 
  #   data: data vector believed to originate from a mixture of unimodal betas.
  #   groups: number of groups/components. 
  #   max_iterations: maximum number of algorithm iterations.
  #   random_initialisation_runs: the number of random runs to be passed to function initialise_me.
  # Returns: 
  #   The number of iterations until convergence.
  #   The parameters (theta) and mixing proportions at final iteration or convergence. 
  #   AIC and BIC for the log likelihood at final iteration or convergence.

  N <- length(data)
  
  theta <- array(0, c(groups, 2, max_iterations + 1))
  mixing_proportions  <- array(0, c(groups, max_iterations + 1))
  log_likelihood <- NULL
  by_group <- NULL
  z <- NULL
  
  initial_values <- initialise_me(data = data, random_runs = random_initialisation_runs)
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
        data = data, # extra argument passed to max_me
        z = z # extra argument passed to max_me
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
    if(iterations == max_iterations + 1) {
      warning("Max iteration limit has been reached without convergence")
    }
  }
  
  AIC <- - 2 * log_likelihood[iterations] + 2 * (3 * groups - 1)
  BIC <- - 2 * log_likelihood[iterations] + (3 * groups - 1) * log(N)   
  
  output <-
    list(
      "Iterations" = iterations,
      "Theta" = theta[, , iterations],
      "MixingProportions" = mixing_proportions[, iterations],
      "AIC" = AIC, 
      "BIC" = BIC,
      "LogLikelihood" = log_likelihood,
      "ClassProbabilities" = z
    )
  class(output) <- "BetaEM"
  return(output)
}

# use microbenchmark to benchmark the function 

## library(microbenchmark)
## microbenchmark(initialise_me(simulated_data), times = 50)
## 
## Unit: seconds
## expr                          min       lq     mean     median uq       max         neval
## run_em(data = simulated_data) 7.36538 8.398256 9.578369 9.1256 10.37926 12.96472    50
