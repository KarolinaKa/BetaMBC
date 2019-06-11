BetaMBC
-----

Model based clustering.
Finite mixtures of unimodal beta distributions (univariate for now). 
EM Algorithm. 

Usage
-----

Simulate the desired mixture in `simulation.R`. 
Follow the example below to calculate density estimates and the clustering solution. 

```r
# Example ///

source("run_em_S3.R")
x <- run_em(simulation$simulated_data)

summary(x) 
plot(x)
beta_clus(x, simulation$components)

```

