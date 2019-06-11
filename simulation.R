##############################################
## Simulation script                        ##
## author: karolina.kulec@ucdconnect.ie     ##
## date: 11 June 2019                       ##
##############################################

source("auxillary_functions.R")

simulation <- simulate_beta(
  N = 1000,
  mixing_proportions = c(0.4, 0.3, 0.3),
  location = c(0, 0.5, 1),
  scale = c(0.03, 0.04, 0.02)
)

 