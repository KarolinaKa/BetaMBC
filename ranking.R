##############################################
## Ranking functions for MultiBetaClus      ##
## author: karolina.kulec@ucdconnect.ie     ##
## date: June 2019                          ##
##############################################

ranking <- function(x,
                    y,
                    RowNames = NULL) 
  {
  # Given two clustering solutions from multi_beta_clus finds observations which belong to different clusters in each.
  # Ranks each observation based on increasing order of clustering uncertainty. 
  # 
  # Args:
  #   x: object of class MultiBetaClus.
  #   y: object of class MultiBetaClus.
  #   RowNames: row names of the original data matrix if available (CpG site identifiers in this case).
  # Returns:
  #   Matrix of class 'ClusRank' with entries ranked in increasing order of uncertainty and possibly uncertainty plots. 
  
  stopifnot(inherits(x, 'MultiBetaClus'))
  stopifnot(inherits(y, 'MultiBetaClus'))
  
  max_uncertainty <- 1 - (1 / x$Groups)
  N <- nrow(x$Data)
  
  uncertainty_matrix1 <-
    cbind(x$ClusterUncertainties, y$ClusterUncertainties)
  membership_matrix <-
    cbind(x$ClusterMemberships, y$ClusterMemberships)
  
  if (!is.null(RowNames)) {
    rownames(uncertainty_matrix1) <- RowNames
    rownames(membership_matrix) <- RowNames
  }
  
  not_identical <-
    which(sapply(1:N, function(i)
      ! identical(membership_matrix[i, 1], membership_matrix[i, 2])))
  
  uncertainty_matrix2 <-
    uncertainty_matrix1[not_identical, ]
  
  output <- uncertainty_matrix2[order(rowSums(uncertainty_matrix2)), ]
  
  class(output) <- "ClusRank"
  return(output)
}

plot.ClusRank <- function(x, Rank, data1, data2) {
  # Plotting function for object of class 'ClusRank'
  #
  # Args:
  #   x: matrix of class ClusRank.
  #   Rank: the desired rank of the CpG site to be plotted.
  #   data1: raw data matrix used to obtain first clustering solution.
  #   data2: raw data matrix used to obtain second clustering solution.
  # Returns:
  #  A histogram of the methylation pattern of the two cohorts at CpG site specified by 'Rank'.
  cohortA <- data1[row.names(x)[Rank], ]
  cohortB <- data2[row.names(x)[Rank], ]
  
  hist(
    as.numeric(cohortA),
    xlim = c(0, 1),
    ylim = c(0, 6),
    freq = F,
    angle = 45,
    density = 90,
    col = "royalblue",
    main = paste("Methylation pattern at CpG site", row.names(cohortA)),
    xlab = expression(beta ~ value)
  )
  hist(
    as.numeric(cohortB),
    add = T,
    freq = F,
    angle = 45,
    density = 80,
    col = "orangered"
  )
  legend(
    "topleft",
    bty = "n",
    legend = c(paste("Cohort A, uncertainty = ", round(as.numeric(x[Rank, 1]), 4)), 
               paste("Cohort B, uncertainty = ", round(as.numeric(x[Rank, 2]), 4))),
    col = c("royalblue", "orangered"),
    pch = 19, 
    cex = 0.9
  )
  
}



