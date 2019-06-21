##############################################
## Simulation script                        ##
## author: karolina.kulec@ucdconnect.ie     ##
## date: 11 June 2019                       ##
##############################################

source("auxillary_functions.R")

# Define location, scale and mixing_proportions vectors each of length 'groups'.

# location <- c(0.2, 0.6, 0.9)
# scale <- c(0.04, 0.08, 0.04)
# mixing_proportions <- c(0.25, 0.5, 0.25)

simulation <- simulate_beta(
  N = 500,
  location,
  scale,
  mixing_proportions,
  groups = 3
)




 