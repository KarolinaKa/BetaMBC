###################################################
## Preliminary analysis                          ##
##                                               ##
## date: March 2019                              ##
## Dataset file: "analysis_data.rds"             ##
###################################################

library(mclust)

analysis_data <- readRDS("analysis_data.rds")

set.seed(420)
analysis_data500 <- 
  analysis_data[sample(nrow(analysis_data), 500),]

mConvert <- function(x)
  # Converts beta values to M values
{
  i <- x / (1 - x)
  log(i, base = 2)
}

betaConvert <- function(x)
  # Converts M values back to beta values
{
  (2 ^ x) / (2 ^ x + 1)
}

analysis_data500 <- mConvert(analysis_data500)

# fit the default Mclust model
data_fit <- Mclust(analysis_data500)
summary(data_fit)
data_fit$BIC

# data_fit plot
par(mfrow = c(3, 3), family = "serif")
sapply(1:data_fit$G, function(i)
  hist(
    betaConvert(data_fit$parameters$mean[, i]),
    xlim = c(0, 1),
    ylab = paste("Cluster", i),
    xlab = "Beta values",
    main = "",
    font.lab = 3
  ))


