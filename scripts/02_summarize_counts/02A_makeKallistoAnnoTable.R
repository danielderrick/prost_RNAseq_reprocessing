# MDD RNASeq Make Kallisto Annotation Table
#
# The purpose of this script is to create an annotation table
# to map Kallisto identifiers (from GENCODE) to ENSEMBL ID.
#
# This script uses sleuth to build the table, but elsewhere DESeq2 is used
# in place of sleuth due to missing features in Sleuth (the ability to
# specify contrasts, etc.)
#
### At some point, this script should be simplified, as it is currently a
# holdover from the initial analysis attempt using sleuth

#####################################################################
# Setup
library(stringr)
library(tidyverse)
library(sleuth)

# Files containing abundance HDF5 files
sample_id <- dir(file.path("~/data/korkola/prostate/", "results"))
kal_dirs <- file.path("~/data/korkola/prostate/", "results", sample_id)
tFixName <- function(x) {
  # One batch has "5nM" in sample names, one does not.
  # This function removes those labels to make names consistent.
  annoError <- grep("nM", x)
  if (exists("annoError")) {
    if (length(annoError) != 0) {
      x <- x[-annoError]
    } else {
      x <- x
    }
  }
  return(x)
}         # To fix inconsistencies in file names

#####################################################################

# Creating metadata file
metadata <- list()
metadata[[1]] <- str_split_fixed(sample_id[1:16], "_", 6)
metadata[[1]] <- metadata[[1]][, -2]

metadata[[2]] <- str_split_fixed(sample_id[17:33], "_", 7)
metadata[[2]] <- metadata[[2]][, -c(2, 4)]

metadata[[3]] <- str_split_fixed(sample_id[34:126], "_", 5)

metadata <- data.frame(Reduce(rbind, metadata))

colnames(metadata) <- c("run", "cell_line", "barcode", "lane", "replicate")
metadata <- metadata %>% dplyr::select(cell_line, lane, replicate)
metadata <- cbind(sample_id, metadata)

# Creating df for sleuth
s2c <- mutate(metadata, path = kal_dirs)
colnames(s2c)[1] <- "sample"

# Slow step - read in files etc.
so <- sleuth_prep(s2c)

so <- sleuth_fit(so, ~cell_line, 'full') # full model = lane + ligand/time
so <- sleuth_fit(so, ~1, 'reduced')       # reduced model = lane
so <- sleuth_lrt(so, 'reduced', 'full')      # compare full model to lane-only

# Getting results table
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
gencode_to_hugo <- data.frame(og_name = sleuth_table$target_id)
temp <- str_split_fixed(gencode_to_hugo$og_name, "[|]", 2)

gencode_to_hugo$txid <- temp[, 1]

# sleuth_table$target_id <- str_extract(sleuth_table$target_id, "ENST[:alnum:]+.[:alnum:]+")
gencode_to_hugo$txid <- str_extract(gencode_to_hugo$txid, "[:alnum:]+")

# Adding HUGO IDs to results table
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'www.ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

colnames(t2g)[1] <- "txid"
t2g2 <- dplyr::inner_join(gencode_to_hugo, t2g, by = "txid")
annoTable <- t2g2
colnames(annoTable)[1] <- "target_id"

dir <- "processed_data/00_prep/"
if (!dir.exists(dir)) {
  dir.create(dir, recursive = TRUE)
}

save(annoTable,
     file = "processed_data/00_prep/kal_annoTable.Rdata")
