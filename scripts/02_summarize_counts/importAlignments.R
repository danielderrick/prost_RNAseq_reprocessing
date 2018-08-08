# MDD RNASeq Abundance to DESeq

# The purpose of this script is to produce gene-level summaries and other
# data for downstream analysis from Kallisto transcript abundances.

###############################################################################
# Setup
library(org.Hs.eg.db)
library(tidyverse)
library(tximport)
library(DESeq2)
library(rhdf5)
library(biomaRt)

# Each sample should have its own directory in dir/ containing an abundance.h5
dir.data <- "data/kal_out/"
load("output/kal_annoTable.Rdata")

###############################################################################
# Creating metadata from file names
files <- list.files(dir.data)
files <- files[-grep("CAL", files)] # removing a breast cancer sample
metadata <- list()
metadata[[1]] <- str_split_fixed(files[1:16], "_", 6)
metadata[[1]] <- metadata[[1]][, -2]

metadata[[2]] <- str_split_fixed(files[17:33], "_", 7)
metadata[[2]] <- metadata[[2]][, -c(2, 4)]

metadata[[3]] <- str_split_fixed(files[34:122], "_", 5)
metadata <- data.frame(Reduce(rbind, metadata))
colnames(metadata) <- c("run", "cell_line", "barcode", "lane", "replicate")
metadata <- dplyr::select(metadata, cell_line, lane, replicate)
meta <- cbind(files, metadata)
meta$cell_line <- as.character(meta$cell_line)
meta[meta$cell_line == "LnCaP", 2] <- "LnCaP-FGC"

# Formatting table to map transcript IDs to gene IDs
txi <- tximport(sprintf("%s%s/%s", dir.data, files, "abundance.h5"),
                type = "kallisto", txOut = TRUE)

# Creating table for mapping to gene-level
tx2g <- 
  data.frame(TXNAME = rownames(txi$abundance)) %>% 
  mutate(tx_id = str_extract(TXNAME, "ENST[:alnum:]+")) %>% 
  dplyr::select(TXNAME, tx_id)

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'uswest.ensembl.org')

tx2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", 
                                      "ensembl_gene_id",
                                      "entrezgene",
                                      "hgnc_symbol"), 
                       mart = mart) %>% 
  dplyr::rename(tx_id = ensembl_transcript_id,
                ensembl_gene = ensembl_gene_id,
                entrez_gene = entrezgene) %>% 
  inner_join(tx2g, by = "tx_id") %>% 
  dplyr::select(TXNAME, ensembl_gene, everything())

# summarizing to gene-level
txi <- summarizeToGene(txi, tx2gene = tx2g[, 1:2])

# Making DESeqDataSet for purpose of extracting counts
meta$cell_line[grep("-", meta$cell_line)] <- str_replace(meta$cell_line[grep("-", meta$cell_line)], "-", "_")
dds.fpkm <- DESeqDataSetFromTximport(txi, meta, ~ 1)

dds.fpkm.collapsed <- collapseReplicates(dds.fpkm, 
                                         groupby = meta$cell_line)
filtro <- rowSums( counts(dds.fpkm.collapsed) >= 5) >= 2
dds.fpkm.collapsed <- dds.fpkm.collapsed[filtro, ]
dds.fpkm.collapsed <- DESeq(dds.fpkm.collapsed, parallel = TRUE)

# Getting FPKM counts and log2 FPKM counts
fpkm.prostate      <- fpkm(dds.fpkm.collapsed)
log2.fpkm.prostate <- apply(fpkm.prostate, c(1,2), function(x) {x <- log2(x + 1)})

# Saving
dir.out <- "output/counts/"
if (!dir.exists(dir.out)) {
  dir.create(dir.out, recursive = TRUE)
}

save(meta, dds.fpkm, fpkm.prostate, log2.fpkm.prostate, 
     file = sprintf("%s%s", dir.out, "prostate_rnaseq_dds_fpkm.Rdata"))

# writing to TSV files

raw.prostate.toWrite <-
  data.frame(counts(dds.fpkm.collapsed)) %>% 
  mutate(ensembl_id = rownames(counts(dds.fpkm.collapsed))) %>% 
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

write.table(raw.prostate.toWrite, 
            file = sprintf("%s%s%s", "output/", "counts/", "prost_raw_counts.tsv"),
            quote = FALSE,
            row.names = FALSE, 
            sep = "\t")

write.table(fpkm.prostate.toWrite, 
            file = sprintf("%s%s%s", "output/", "counts/", "prost_fpkm.tsv"),
            quote = FALSE,
            row.names = FALSE, 
            sep = "\t")

write.table(log2.fpkm.prostate.toWrite, 
            file = sprintf("%s%s%s", "output/", "counts/", "prost_log2fpkm.tsv"),
            quote = FALSE,
            row.names = FALSE, 
            sep = "\t")


