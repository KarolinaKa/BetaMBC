##############################################
## Auxillary functions for the EM algorithm ##
##                                          ##  
##  1. simulate_beta                        ##  
##  2. dbeta_rep                            ##
##  3. initialise_me                        ##
##  4. location_constraint and its inverse  ##
##  5. max_me                               ##  
##                                          ##
## Other functions                          ##  
##                                          ##
##  1. adjusted_rand                        ##
##  3. permute_me                           ##
##                                          ##
## author: karolina.kulec@ucdconnect.ie     ##
## date: June 2019                          ##
##############################################

simulate_beta <-
  function(N,
           mixing_proportions,
           location,
           scale,
           groups = 3) {
  # Returns random draws from the reparametrised beta density.
  #
  # Args: 
  #   N: number of draws required.
  #   mixing_proportions: mixing proportions vector of length groups.
  #   location: mode vector of length groups.
  #   scale: spread parameter (analagous to variance in Gaussian setting) vector of length groups.
  #   groups: number of groups.
  # Returns: 
  #   Vector of components.
  #   Random draws from the reparametrised beta density.
    
    # error handling
    if (length(mixing_proportions) != groups) stop(paste("mixing_proportions must have length", groups))
    if (length(location) != groups) stop(paste("location must have length", groups))
    if (length(scale) != groups) stop(paste("scale must have length", groups))
    
    components <-
      sample(1:groups,
             prob = mixing_proportions,
             size = N,
             replace = TRUE)
    simulated_data <-
      rbeta(N, location[components] / scale[components] + 1, (1 - location[components]) / scale[components] + 1)
    list(components = components, simulated_data = simulated_data)
  }

dbeta.rep <-
  function(x,
           location,
           scale) {
    # Defines the density of the reparametrised beta function using the inbuilt R dbeta function.
    #
    # Args: 
    #   x: evaluation points in [0, 1].
    #   location: mode. 
    #   scale: spread parameter (analagous to variance in Gaussian setting).
    # Returns: 
    #   Reparametrised beta distribution density.
    
    dbeta(x, location / scale + 1, (1 - location) / scale + 1)
  }

initialise_me <- function(data,
                          random_runs,
                          groups = 3,
                          scale_min = 0.01,
                          scale_max = 0.9) {
  # Random initialisation function for the EM algorithm.
  #
  # Args: 
  #   data: data vector generated using a finitie mixture of unimodal beta distributions.
  #   random_runs: number of random initialisations.
  #   groups: number of groups.
  #   scale_min and scale_max: min and max of the interval from which to sample the scale parameter.
  # Returns: 
  #   A list containing two elements: (1) a parameter matrix and (2) a mixing proportions vector.
  
  N <- length(data)
  
  theta <- array(0, c(groups, 2, random_runs))
  final_theta <- array(0, c(groups, 2, 1))
  
  random_numbers <- array(0, c(groups, 1, random_runs))
  mixing_proportions <- array(0, c(groups, 1, random_runs))
  final_mixing_proportions <- array(0, c(groups, 1, 1))
  
  log_likelihood <- NULL
  
  for (i in 1:random_runs) {
    random_numbers[, , i] <- runif(groups, 0, 1)
    mixing_proportions[, , i] <-
      random_numbers[, , i] / sum(random_numbers[, , i])
    theta[, 1, i] <- kmeans(data, groups)$centers
    theta[, 2, i] <- rep(runif(1, scale_min, scale_max), groups)
    
    log_likelihood[i] <- sum(log(
      sapply(1:groups, function(j)
        dbeta.rep(data, theta[j, 1, i], theta[j, 2, i])) %*% mixing_proportions[, , i]
    ))
  }
  
  index <- which(log_likelihood == max(log_likelihood))
  final_mixing_proportions <- mixing_proportions[, , index]
  final_theta <- theta[, , index]
  
  output <- list(final_theta, final_mixing_proportions)
}

# use microbenchmark to benchmark the function 

## library(microbenchmark)
## microbenchmark(initialise_me(simulated_data), times = 50)
##
## Unit: milliseconds
## expr                          min      lq       mean     median   uq       max      neval
## initialise_me(simulated_data) 299.5482 324.9853 349.5559 338.8011 367.4486 561.5888 50

location_constraint <- function(x) {
  # Used during optimisation ensures location parameters in [0, 1]. Similar to sigmoid curve.
  # 
  # Args: 
  #   x: evaluation points.
  # Returns: 
  #   Values in [0, 1].
  
  (1 / pi) * atan(x) + 1 / 2
}

inv_location_constraint <- function(x) {
  # Inverse of location_constraint.
  
  tan((x - 1 / 2) * pi)
}

max_me <- function(theta, data, z, groups = 3) {
  # Defines the objective function to be maximised.
  #
  # Args: 
  #   theta: vector of parameters to be maximised over.
  #   data: data vector generated using a finitie mixture of unimodal beta distributions. 
  #   z: (N x groups) posterior probability matrix.
  #   groups: number of components/groups in the mixture. 
  # Returns: 
  #   Negative objective function output (because optim minimises by default).
  
  ind <-
    sum(sapply(1:groups, function(j)
      z[, j] %*% log(
        dbeta.rep(
          data,
          location = location_constraint(theta[2 * j - 1]),
          scale = exp(theta[2 * j])
        )
      )))
  return(-ind)                            
}

adjusted_rand <- function(table) {
  # Calculates the Adjusted Rand Index (Wikipedia definition) for two clustering solutions.
  #
  # Args:
  #   table: cross tabulation of two clustering solutions.
  # Returns:
  #   Adjusted Rand Index.
  N <- sum(table)
  n <- diag(table)
  a <- rowSums(table)
  b <- colSums(table)
  
  numerator <-
    sum(choose(n, 2)) - (sum(choose(a, 2)) *  sum(choose(b, 2))) / choose(N, 2)
  denominator <-
    1 / 2 * (sum(choose(a, 2)) + sum(choose(b, 2))) - (sum(choose(a, 2)) * sum(choose(b, 2))) / choose(N, 2)
  
  numerator / denominator
}

permute_me <- function(x) {
  # Calculates all vector permutations.
  #
  # Args:
  #   x: number of groups used to create the vector 1:x.
  # Returns:
  #   Permutation matrix (number of permutations x groups).
  combination_matrix <-
    expand.grid(1:x, 1:x, 1:x)
  permutation_matrix <-
    combination_matrix[apply(combination_matrix, 1, function(i)
      length(unique(i)) == 3),]
  
  row.names(permutation_matrix) <- NULL
  colnames(permutation_matrix) <- NULL
  
  as.matrix(permutation_matrix)
}

