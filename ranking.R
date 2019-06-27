##############################################
## Ranking functions for MultiBetaClus      ##
## author: karolina.kulec@ucdconnect.ie     ##
## date: June 2019                          ##
##############################################

# Code is semi-specific to the DNA methylation data provided in the repo. 

ranking <- function(x,
                    y,
                    RowNames = NULL,
                    plot = TRUE) {
  # Given two clustering solutions from multi_beta_clus finds observations which belong to different clusters in each.
  # Ranks each observation based in increasing order of clustering uncertainty. 
  # 
  # Args:
  #   x: object of class MultiBetaClus.
  #   y: object of class MultiBetaClus.
  #   RowNames: row names of the original data matrix if available (CpG site identifiers in this case).
  #   plot: specify if plots are required. 
  # Returns:
  #   Matrix of class 'ClusRank' with entries ranked in increasing order of uncertainty.
  #   Uncertainty plots if plot = TRUE.
  
  stopifnot(inherits(x, 'MultiBetaClus'))
  stopifnot(inherits(y, 'MultiBetaClus'))
  
  max_uncertainty <- 1 - (1 / x$Groups)
  N <- nrow(x$Data)
  
  uncertainty_matrix <-
    cbind(x$ClusterUncertainties, y$ClusterUncertainties)
  membership_matrix <-
    cbind(x$ClusterMemberships, y$ClusterMemberships)
  
  if (!is.null(RowNames)) {
    rownames(uncertainty_matrix) <- RowNames
    rownames(membership_matrix) <- RowNames
  }
  
  different_cpg <-
    which(sapply(1:N, function(i)
      ! identical(membership_matrix[i, 1], membership_matrix[i, 2])))
  
  different_cpg_uncertainties <-
    uncertainty_matrix[different_cpg, ]
  sort_cpg_uncertainties <-
    different_cpg_uncertainties[order(rowSums(different_cpg_uncertainties)), ]
  
  if (plot == TRUE) {
    n <- length(different_cpg)
    par(mfrow = c(1, 2))
    plot(
      1:n,
      sort(sort_cpg_uncertainties[, 1]),
      type = "h",
      ylim = c(0, 1),
      xlab = "CpG site",
      ylab = "Uncertainty",
      main = "Basal-like Cohort"
    )
    abline(h = max_uncertainty, col = "red")
    text(100, max_uncertainty + 0.03, labels = c("Maximum Uncertainty"), cex = 0.8)
    grid()
    
    plot(
      1:n,
      sort(sort_cpg_uncertainties[, 2]),
      type = "h",
      ylim = c(0, 1),
      xlab = "CpG site",
      ylab = "Uncertainty",
      main = "Non Basal-like Cohort"
    )
    abline(h = max_uncertainty, col = "red")
    text(100, max_uncertainty + 0.03, labels = c("Maximum Uncertainty"), cex = 0.8)
    grid()
    
    par(mfrow = c(1, 1))
    plot(1:n, sort(abs(different_cpg_uncertainties[, 1] - different_cpg_uncertainties[, 2])), type = "h", ylim = c(0, 1),
         xlab = "CpG site",
         ylab = "Uncertainty",
         main = "Absolute difference in uncertainty between two cohorts")
    abline(h = max_uncertainty, col = "red")
    text(100, max_uncertainty + 0.03, labels = c("Maximum Uncertainty"), cex = 0.8)
    grid()
  }
  class(sort_cpg_uncertainties) <- "ClusRank"
  return(sort_cpg_uncertainties)
}


plot.ClusRank <- function(x, Rank, data1, data2) {
  # Plotting function for 'ranking' 
  #
  # Args:
  #   x: matrix of class ClusRank.
  #   Rank: the desired rank of the CpG site to be plotted.
  #   data1: raw data matrix used to obtain first clustering solution.
  #   data2: raw data matrix used to obtain second clustering solution.
  # Returns:
  #  A histogram of the methylation pattern of the two cohorts at CpG site specified by 'Rank'.
  basal_like <- data1[row.names(x)[Rank], ]
  non_basal_like <- data2[row.names(x)[Rank], ]
  
  hist(
    as.numeric(basal_like),
    xlim = c(0, 1),
    freq = F,
    angle = 45,
    density = 90,
    breaks =  25, 
    col = "royalblue",
    main = paste("Methylation pattern at CpG site", row.names(basal_like)),
    xlab = expression(beta ~ value)
  )
  hist(
    as.numeric(non_basal_like),
    add = T,
    freq = F,
    breaks =  25, 
    angle = 45,
    density = 80,
    col = "orangered"
  )
  legend(
    "topright",
    bty = "n",
    legend = c(paste("Basal-Like, uncertainty = ", round(as.numeric(x[Rank, 1]), 4)), 
               paste("Non Basal-Like, uncertainty = ", round(as.numeric(x[Rank, 2]), 4))),
    col = c("royalblue", "orangered"),
    pch = 19, 
    cex = 0.9
  )
  
}


