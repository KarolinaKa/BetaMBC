###############################################
## TCGA-BRCA Data prep                       ##
## date: June 2019                           ##
## Dataset file: "BRCA.methylation.27k.450k" ##
## Supplementary file: "Table1Nature.csv"    ##
###############################################

# Load the DNA methylation dataset
methylation_data    <-
  read.table("BRCA.methylation.27k.450k.txt", header = TRUE)
methylation_data_ID <- substr(names(methylation_data), 1, 12)
methylation_data    <- as.matrix(methylation_data)

# Find and remove rows with missing values
missing_methylation_data <-
  which(is.na(methylation_data), arr.ind = T)
uniquely_missing         <-
  unique(row.names(missing_methylation_data))
methylation_data         <-
  methylation_data[!row.names(methylation_data) %in% uniquely_missing,]

# Load the reference table which relates subject ID to clinical profile
reference_table <- read.csv("Table1Nature.csv", header = TRUE)

# Format matching ID's
reference_ID        <- reference_table$Complete.TCGA.ID
reference_ID        <- gsub("-", ".", as.character(reference_ID))

# Math the subject with subtype
match_us           <-
  match(methylation_data_ID, reference_ID, nomatch = NA)
reference_subtypes <- reference_table$PAM50.mRNA
reference_subtypes <- reference_subtypes[match_us]

methylation_data   <-
  methylation_data[, !is.na(reference_subtypes)]
reference_subtypes <-
  reference_subtypes[!is.na(reference_subtypes)]

# save the files in compressed format
saveRDS(methylation_data, "analysis_data.rds")
saveRDS(reference_subtypes, "analysis_subtypes.rds")
