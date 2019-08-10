##############################################
## Ranking functions for MultiBetaClus      ##
## author: karolina.kulec@ucdconnect.ie     ##
## date: June 2019                          ##
##############################################

ranking <- function(x,
                    y,
                    RowNames = NULL) 
  {
  # Given two clustering solutions from multi_beta_clus 
  # finds observations which belong to different clusters in each
  # and ranks each observation based on increasing order of 
  # clustering uncertainty. 
  # 
  # Args:
  #   x: object of class MultiBetaClus.
  #   y: object of class MultiBetaClus.
  #   RowNames: row names of the original data matrix 
  #   if available (CpG site identifiers in this case).
  # Returns:
  #   Matrix of class 'ClusRank' with entries ranked in increasing 
  #   order of uncertainty. 
  
  stopifnot(inherits(x, 'MultiBetaClus'))
  stopifnot(inherits(y, 'MultiBetaClus'))
  
  N <- nrow(x$Data)
  
  uncertainty_matrix <-
    cbind(x$ClusterUncertainties, y$ClusterUncertainties)
  membership_matrix <-
    cbind(x$ClusterMemberships, y$ClusterMemberships)
  
  if (!is.null(RowNames)) {
    rownames(uncertainty_matrix) <- RowNames
    rownames(membership_matrix) <- RowNames
  }
  
  not_identical <-
    which(sapply(1:N, function(i)
       ! identical(membership_matrix[i, 1], 
                  membership_matrix[i, 2])))
  
  uncertainty_matrix <-
    uncertainty_matrix[not_identical, ]
  
  output <- 
    uncertainty_matrix[order(rowSums(uncertainty_matrix)), ]
  
  class(output) <- "ClusRank"
  return(output)
}

plot.ClusRank <- function(x, Rank, data1, data2, ylim) {
  # Plotting function for object of class 'ClusRank'
  #
  # Args:
  #   x: matrix of class ClusRank.
  #   Rank: the desired rank of the CpG site to be plotted.
  #   data1: raw data matrix used to obtain first clustering solution.
  #   data2: raw data matrix used to obtain second clustering solution.
  #   ylim: sets the upper ylim for plots
  # Returns:
  #  A histogram of the methylation pattern of the two cohorts at CpG site specified by 'Rank'.
  cohortA <- data1[row.names(x)[Rank], ]
  cohortB <- data2[row.names(x)[Rank], ]
  
  hist(
    as.numeric(cohortA),
    xlim = c(0, 1),
    ylim = c(0, ylim),
    freq = F,
    angle = 45,
    density = 70,
    col = "royalblue",
    main = paste(row.names(x)[Rank], "rank", Rank),
    xlab = "Beta value", 
    font.lab = 3, 
    font.main = 4, 
    ylab = "Density"
  ) 
  hist(
    as.numeric(cohortB),
    add = T,
    freq = F,
    angle = 45,
    density = 70,
    col = "orangered"
  )
  
}



