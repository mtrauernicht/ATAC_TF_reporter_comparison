# Libraries ---------------------------------------------------------------

library(here)
library(glue)
library(data.table)
library(readr)

library(GenomicRanges)
library(rtracklayer)
library(IRanges)
library(GenomeInfoDb)
library(dplyr)


# Directories and data ----------------------------------------------------

files_df <- readRDS("/DATA/usr/m.trauernicht/projects/ATAC_TF_reporter_comparison/ATAC_seq_analyses/rds/bamfile_atacseq_metadata2_human.rds")
bw_dir <- "/DATA/usr/m.trauernicht/projects/ATAC_TF_reporter_comparison/ATAC_seq_analyses/bigwig/20240626_human"
size_factors <- read_csv("/DATA/usr/m.trauernicht/projects/ATAC_TF_reporter_comparison/ATAC_seq_analyses/bigwig/size_factors_deseq.csv") %>%
  mutate(sample = paste(celltype, treatment, replicate, sep = "_"))

# Split up pools ----------------------------------------------------------


files_df <- files_df[files_df$run == "7822",]

files <- split(files_df$tabix_file, interaction(files_df$celltype, files_df$treatment, files_df$replicate, sep = "_"))
files <- files[lapply(files,length)>0]

seqinfo <- SeqinfoForUCSCGenome("hg38")
seqinfo <- keepStandardChromosomes(seqinfo, "Homo_sapiens")
seqinfo <- seqinfo[names(seqinfo)[-length(seqinfo)]]
seqnames <- seqnames(seqinfo)
gsize <- sum(seqlengths(seqinfo))

for (group in names(files)) {
  f <- files[[group]]
  outfile <- paste0(bw_dir, group, "_3.bw")
  if (file.exists(outfile)) {
    next
  }
  
  data <- lapply(f, fread)
  data <- rbindlist(data)

  
  data <- with(data, GRanges(V1, IRanges(V2, V3)))
  
  factor <- size_factors$sizefactor[size_factors$sample==group]
  
  data <- coverage(data) / factor
  data <- data[names(data) %in% seqnames]
  export.bw(data, paste0("/DATA/usr/m.trauernicht/projects/ATAC_TF_reporter_comparison/ATAC_seq_analyses/bigwig/", group, "_3.bw"))
}
