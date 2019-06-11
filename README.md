# MScThesis2019

Model based clustering.
Finite mixtures of unimodal beta distributions. 
EM Algorithm. 

Usage
-----

Simulate the desired mixture in `simulation.R`. 

```r
# Example ///

source("run_em_S3.R")
x <- run_em(data = simulation$simulated_data)

summary(x) 
plot(x)
beta_clus(x)

```

