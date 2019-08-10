#################################################
## Basal-like versus remaining subtypes        ##
## Luminal A versus Luminal B                  ##
##                                             ##  
## date: July 2019                             ##
## Dataset file: "analysis_data.rds"           ##
## Supplementary file: "analysis_subtypes.rds" ##
#################################################

source("multi_beta_clus.R")
source("ranking.R")

analysis_data <- readRDS("analysis_data.rds")
analysis_subtypes <- readRDS("analysis_subtypes.rds")

basal_like <-
  analysis_data[, analysis_subtypes == "Basal-like"]
non_basal_like <-
  analysis_data[, analysis_subtypes != "Basal-like"]

luminal_A <- 
  analysis_data[, analysis_subtypes == "Luminal A"]
luminal_B <- 
  analysis_data[, analysis_subtypes == "Luminal B"]

################################
## Luminal A versus Luminal B ##
################################

mbc_luma <- multi_beta_clus(
  luminal_A,
  groups = 3,
  location = c(0, 0.5, 1),
  scale = c(0.01, 0.1, 0.1)
)

mbc_lumb <- multi_beta_clus(
  luminal_B,
  groups = 3,
  location = c(0, 0.5, 1),
  scale = c(0.01, 0.1, 0.1)
)

table(mbc_luma$ClusterMemberships)
table(mbc_lumb$ClusterMemberships)

tab <- table(
  "Luminal A" = mbc_luma$ClusterMemberships,
  "Luminal B" = mbc_lumb$ClusterMemberships
)

cluster <- 
  cbind("Luminal A" = mbc_luma$ClusterMemberships, 
        "Luminal B" = mbc_lumb$ClusterMemberships)
which(cluster[, "Luminal A"] == 1 &
        cluster[, "Luminal B"] == 3)

# Cluster uncertainties 
mbc_luma$ClusterUncertainties[2609]
mbc_lumb$ClusterUncertainties[2609]

mbc_luma$ClusterUncertainties[17206]
mbc_lumb$ClusterUncertainties[17206]

par(mfrow = c(1, 2), 
    family = "serif", 
    font.lab = 3, 
    font.main = 4)

hist(
  luminal_A[2609,],
  xlim = c(0, 1),
  density = 80,
  angle = 65,
  freq = F,
  col = "royalblue",
  xlab = "Beta value",
  ylim = c(0, 4),
  main = "cg03544320"
)
hist(
  luminal_B[2609,],
  add = T,
  density = 80,
  angle = 65,
  freq = F,
  col = "orangered"
)

hist(
  luminal_A[17206,],
  xlim = c(0, 1),
  density = 80,
  angle = 65,
  freq = F,
  col = "royalblue",
  xlab = "Beta value",
  ylim = c(0, 4),
  main = "cg23391785"
)
hist(
  luminal_B[17206,],
  add = T,
  density = 80,
  angle = 65,
  freq = F,
  col = "orangered"
)

##########################################
## Basal-like versus remaining subtypes ##
##########################################

mbc_bl <- multi_beta_clus(
  basal_like,
  groups = 3,
  location = c(0, 0.5, 1),
  scale = c(0.01, 0.1, 0.1)
)

mbc_nbl <- multi_beta_clus(
  non_basal_like,
  groups = 3,
  location = c(0, 0.5, 1),
  scale = c(0.01, 0.1, 0.1)
)

table(mbc_nbl$ClusterMemberships)

saveRDS(mbc_bl, "mbc_bl.rds")
saveRDS(mbc_nbl, "mbc_nbl.rds")

# cross-tabulation of results
tab <- table(
  "basal_like" = mbc_bl$ClusterMemberships,
  "non_basal_like" = mbc_nbl$ClusterMemberships
)

ranks <- ranking(mbc_bl, mbc_nbl,
                 RowNames = rownames(basal_like))


######################
## Chi squared test ##
######################

props_lum <- matrix(
  data = c(1529, 804, 18774, 19499), 
  ncol = 2, nrow = 2, byrow = T)
prop.test(props_lum, alternative = "greater")


##################################################
## Comparison with other methods                ##
##                                              ##
## date: July 2019                              ##
## Supplementary file: "TCGA_ME-DE_results.csv" ##
##################################################

#######################
## two-sample t-test ##
#######################

# M value transformation
mConvert <- function(x)
{
  i <- x / (1 - x)
  log(i, base = 2)
}

bl_m <- mConvert(basal_like)
nbl_m <- mConvert(non_basal_like)

# Carry out the t-test
tested <- data.frame("p_value" =
                       sapply(1:nrow(bl_m), function(i)
                         t.test(bl_m[i,],
                                nbl_m[i,])$p.value))
rownames(tested) <- row.names(bl_m)

significance <- 0.001

# Adjust p-values for multiple comparisons
tested$BH <- p.adjust(tested$p_value, method = "BH")
tested <- tested[order(tested$BH), ]
tested <- subset(tested, tested$BH < significance)
dim(tested)

# Results compare
compare_tested <-
  match(row.names(ranks), row.names(tested))
sum(!is.na(compare_tested)) /
  length(compare_tested) # % matches

####################################
## Lock and Dunson (2015) results ##
####################################

ld_results <-
  read.csv("TCGA_ME-GE_results.csv", header = TRUE)
ld_results <-
  ld_results[-which(which(is.na(ld_results),
                          arr.ind = T)[, 2] == 3), ]
ld_results <- ld_results[order(ld_results$ME_Posterior_Prob), ]
ld_results <-
  subset(ld_results, ld_results$ME_Posterior_Prob < significance)
dim(ld_results)

# Results compare
compare_ld <- match(row.names(ranks), ld_results$CpG)
sum(!is.na(compare_ld)) / length(compare_ld) # % matches
dim(ld_results)

