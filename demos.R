source("run_em_S3.R")

# Real Data Demo 

NonBasalLike <- read.csv("TCGA_NonBasalLike.csv", row.names = 1)

# One simple denisty estimation (one and then 2 clusters)
test_data <- as.vector(t(NonBasalLike[13,]))

x <- run_em(test_data, groups = 1)
summary(x)
plot(x)

x2 <- run_em(test_data, groups = 2)
summary(x2)
plot(x2)
clus_solution <- beta_clus(x2)

# Possible 3 clusters 
test_data <- as.vector(t(NonBasalLike[81,]))

x3 <- run_em(test_data, groups = 3, max_iterations = 500, convergence_limit = 0.001)
summary(x3)
plot(x3)
clus_solution <- beta_clus(x3)

# Simulation Demo

## Two Clusters
## Easy 

# simulation <- simulate_beta(
#  N = 500,
#  location = c(0.2, 0.8),
#  scale = c(0.04, 0.08), 
#  mixing_proportions = c(0.45, 0.55),
#  groups = 2
# )

x <- run_em(simulation$SimulatedData, groups = 2)
summary(x)
plot(x)
simulation_plot(x, simulation)
label_switch(x, simulation)
clus_solution <- simulation_beta_clus(x, simulation)
clus_solution$CrossTab
clus_solution$AdjustedRandIndex

## Simple Example 3 clusters
## Perfect Cluster Separation 

# simulation <- simulate_beta(
#  N = 500,
#  location = c(0.15, 0.5, 0.9),
#  scale = c(0.01, 0.02, 0.01), 
#  mixing_proportions = c(0.3, 0.5, 0.2),
#  groups = 3
# )

x <- run_em(simulation$SimulatedData, groups = 3)
summary(x)
plot(x)
simulation_plot(x, simulation)
label_switch(x, simulation)
clus_solution <- simulation_beta_clus(x, simulation)
clus_solution$CrossTab
clus_solution$AdjustedRandIndex

## Tough Example 

# simulation <- simulate_beta(
# N = 500,
# location = c(0.2, 0.5, 0.8),
# scale = c(0.04, 0.08, 0.04), 
# mixing_proportions = c(0.3, 0.5, 0.2),
# groups = 3
# )

x <- run_em(simulation$SimulatedData, groups = 3, max_iterations = 500)
summary(x)
plot(x)
simulation_plot(x, simulation)
label_switch(x, simulation)
clus_solution <- simulation_beta_clus(x, simulation)
clus_solution$CrossTab
clus_solution$AdjustedRandIndex






