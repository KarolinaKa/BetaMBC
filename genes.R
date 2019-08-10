#######################################################
## Genes                                             ##
##                                                   ##
## date: July 2019                                   ##
## Dataset file: "analysis_data.rds"                 ##
## Supplementary files: "analysis_subtypes.rds"      ##
##                      "TCGA_ME-GE_results.csv.rds" ##
##                      "mbc_bl.rds"                 ##
##                      "mbc_nbl.rds"                ##
#######################################################

mbc_bl <- readRDS("mbc_bl.rds")
mbc_nbl <- readRDS("mbc_nbl.rds")

analysis_data <- readRDS("analysis_data.rds")
analysis_subtypes <- readRDS("analysis_subtypes.rds")

basal_like <-
  analysis_data[, analysis_subtypes == "Basal-like"]
non_basal_like <-
  analysis_data[, analysis_subtypes != "Basal-like"]

ld_results <-
  read.csv("TCGA_ME-GE_results.csv", header = TRUE)
ld_results$Genes <- as.character(ld_results$Genes)


# Genes of interest with correponding CpG sites

genes_of_interest <-
  c("MGMT", "MLH1", "DAPK1", "CDH1", "CDH13", "GSTP1", "SFN")
goi <-
  ld_results[which(ld_results$Genes %in% genes_of_interest), 1:2]

# Basal-like
table(mbc_bl$ClusterMemberships[which(row.names(basal_like) %in% goi$CpG)]) 
table(mbc_bl$ClusterMemberships[-which(row.names(basal_like) %in% goi$CpG)]) 

props_bl <- matrix(data = c(14, 36, 2923, 17330), ncol = 2, nrow = 2, byrow = T)
fisher.test(props_bl)

# Non-Basal-like
table(mbc_nbl$ClusterMemberships[which(row.names(non_basal_like) %in% goi$CpG)]) 
table(mbc_nbl$ClusterMemberships[-which(row.names(non_basal_like) %in% goi$CpG)])  

props_nbl <- matrix(data = c(15, 35, 3034, 17219), ncol = 2, nrow = 2, byrow = T)
fisher.test(props_nbl)

##########

bl_goi_stats <-
  data.frame(
    "CpG" = row.names(basal_like)[which(row.names(basal_like) %in% goi$CpG)],
    "Genes" = ld_results[which(as.character(ld_results$CpG) %in% 
                                 row.names(basal_like)[which(row.names(basal_like) %in% 
                                                               goi$CpG)]), 2],
    "ClusterMembership" = mbc_bl$ClusterMemberships[which(row.names(basal_like) %in% goi$CpG)],
    "ClusterUncertainty" = mbc_bl$ClusterUncertainties[which(row.names(basal_like) %in% goi$CpG)]
  )


table(bl_goi_stats[, 2:3])

nbl_goi_stats <-
  data.frame(
    "CpG" = row.names(non_basal_like)[which(row.names(non_basal_like) %in% goi$CpG)],
    "Genes" = ld_results[which(as.character(ld_results$CpG) %in% 
                                 row.names(non_basal_like)[which(row.names(non_basal_like) %in% 
                                                               goi$CpG)]), 2],
    "ClusterMembership" = mbc_nbl$ClusterMemberships[which(row.names(non_basal_like) %in% goi$CpG)],
    "ClusterUncertainty" = mbc_nbl$ClusterUncertainties[which(row.names(non_basal_like) %in% goi$CpG)]
  )

table(nbl_goi_stats[, 2:3])



subset(nbl_goi_stats, nbl_goi_stats$ClusterMembership == 3)

hist(basal_like["cg16777782",], xlim = c(0, 1))
hist(non_basal_like["cg16777782", ], xlim = c(0, 1))

mbc_bl$ClusterMemberships[which(row.names(basal_like) == "cg16777782")]


