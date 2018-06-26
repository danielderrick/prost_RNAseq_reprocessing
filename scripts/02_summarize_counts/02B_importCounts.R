# MDD RNASeq Abundance to DESeq

# The purpose of this script is to produce gene-level summaries and other
# data for downstream analysis from Kallisto transcript abundances.

###############################################################################
# Setup
library(tidyverse)
library(tximport)
library(DESeq2)
library(rhdf5)

# Each sample should have its own directory in dir/ containing an abundance.h5
dir <- "~/data/korkola/prostate/results/"
load("processed_data/00_prep/kal_annoTable.Rdata")
###############################################################################

# Creating metadata from file names
samples <- list.files(dir)

metadata <- list()
metadata[[1]] <- str_split_fixed(samples[1:16], "_", 6)
metadata[[1]] <- metadata[[1]][, -2]

metadata[[2]] <- str_split_fixed(samples[17:33], "_", 7)
metadata[[2]] <- metadata[[2]][, -c(2, 4)]

metadata[[3]] <- str_split_fixed(samples[34:126], "_", 5)
metadata <- data.frame(Reduce(rbind, metadata))
colnames(metadata) <- c("run", "cell_line", "barcode", "lane", "replicate")
metadata <- metadata %>% dplyr::select(cell_line, lane, replicate)
meta <- cbind(samples, metadata)
meta$cell_line <- as.character(meta$cell_line)
meta[meta$cell_line == "LnCaP", 2] <- "LnCaP-FGC"

# Formatting table to map transcript IDs to gene IDs
annoTable <- 
  annoTable %>% 
  dplyr::rename(TXNAME = target_id, GENEID = ens_gene) %>%
  dplyr::select(TXNAME, GENEID)

# excluding the one breast sample that was with these files
files <- file.path(dir, meta$samples, "abundance.h5")
names(files) <- meta$samples
files <- files[-grep("CAL", names(files))]
meta <- meta[-grep("CAL", meta$samples), ]

# Importing files
txi <- tximport(files, type = "kallisto", tx2gene = annoTable)

# Making DESeqDataSet for purpose of extracting counts
dds.fpkm <- DESeqDataSetFromTximport(txi, meta, ~ cell_line)

dds.fpkm.collapsed <- collapseReplicates(dds.fpkm, groupby = meta$cell_line)
dds.fpkm <- DESeq(dds.fpkm.collapsed)

# Getting FPKM counts and log2 FPKM counts
fpkm.prostate      <- fpkm(dds.fpkm)
log2.fpkm.prostate <- apply(fpkm.prostate, c(1,2), function(x) {x <- log2(x + 1)})

# Saving
dir <- "output/"
if (!dir.exists(dir)) {
  dir.create(dir, recursive = TRUE)
}
if (!dir.exists(sprintf("%s%s", dir, "counts/", recursive = TRUE))) {
  dir.create(sprintf("%s%s", dir, "counts/", recursive = TRUE))
}

save(meta, dds.fpkm, fpkm.prostate, log2.fpkm.prostate, 
     file = sprintf("%s%s", dir, "prostate_rnaseq_dds_fpkm.Rdata"))

# writing to TSV files

raw.prostate.toWrite <-
  data.frame(counts(dds.fpkm)) %>% 
  mutate(ensembl_id = rownames(counts(dds.fpkm))) %>% 
  mutate(gene_symbol = mapIds(org.Hs.eg.db, ensembl_id, "SYMBOL", "ENSEMBL")) %>% 
  dplyr::select(ensembl_id, gene_symbol, everything())

fpkm.prostate.toWrite <-
  data.frame(fpkm.prostate) %>% 
  mutate(ensembl_id = rownames(fpkm.prostate)) %>% 
  mutate(gene_symbol = mapIds(org.Hs.eg.db, ensembl_id, "SYMBOL", "ENSEMBL")) %>% 
  dplyr::select(ensembl_id, gene_symbol, everything())

log2.fpkm.prostate.toWrite <-
  data.frame(log2.fpkm.prostate) %>% 
  mutate(ensembl_id = rownames(log2.fpkm.prostate)) %>% 
  mutate(gene_symbol = mapIds(org.Hs.eg.db, ensembl_id, "SYMBOL", "ENSEMBL")) %>% 
  dplyr::select(ensembl_id, gene_symbol, everything())


write.table(fpkm.prostate.toWrite, 
            file = sprintf("%s%s%s", dir, "counts/", "prost_fpkm.tsv"),
            quote = FALSE,
            row.names = FALSE, 
            sep = "\t")

write.table(log2.fpkm.prostate.toWrite, 
            file = sprintf("%s%s%s", dir, "counts/", "prost_log2fpkm.tsv"),
            quote = FALSE,
            row.names = FALSE, 
            sep = "\t")

write.table(raw.prostate.toWrite, 
            file = sprintf("%s%s%s", dir, "counts/", "prost_raw_counts.tsv"),
            quote = FALSE,
            row.names = FALSE, 
            sep = "\t")
