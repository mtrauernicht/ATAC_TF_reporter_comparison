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

files_df <- readRDS("/DATA/usr/m.trauernicht/projects/ATAC_TF_reporter_comparison/ATAC_seq_analyses/rds/bamfile_atacseq_metadata_mES_selected_mt20260420.rds")
bw_dir <- "/DATA/usr/m.trauernicht/projects/ATAC_TF_reporter_comparison/ATAC_seq_analyses/bigwig/mouse_selected_mt20260420"
size_factors <- read_csv("/DATA/usr/m.trauernicht/projects/ATAC_TF_reporter_comparison/ATAC_seq_analyses/bigwig/size_factors_deseq_mt20260420.csv") %>%
  mutate(sample = paste(celltype, treatment, replicate, sep = "_"))

# Split up pools ----------------------------------------------------------


files_df <- files_df[files_df$run %in% c("7784", "7822", "8010", "8011"),]
files_df$replicate[files_df$celltype == "mES" & files_df$treatment == "2i_LIF" & files_df$replicate == "R3" & files_df$run == "8011"] <- "R4"
files_df <- files_df[!(files_df$treatment == "2i_LIF" & files_df$replicate == "R3"), ]

files <- split(files_df$tabix_file, interaction(files_df$celltype, files_df$treatment, files_df$replicate, sep = "_"))
files <- files[lapply(files,length)>0]

seqinfo <- SeqinfoForUCSCGenome("mm10")
seqinfo <- keepStandardChromosomes(seqinfo, "Mus_musculus")
seqinfo <- seqinfo[names(seqinfo)[-length(seqinfo)]]
seqnames <- seqnames(seqinfo)
gsize <- sum(seqlengths(seqinfo))

for (group in names(files)) {
  f <- files[[group]]
  outfile <- paste0(bw_dir, group, ".bw")
  if (file.exists(outfile)) {
    next
  }
  
  data <- lapply(f, fread)
  data <- rbindlist(data)

  
  data <- with(data, GRanges(V1, IRanges(V2, V3)))
  
  factor <- size_factors$sizefactor[size_factors$sample==group]
  
  data <- coverage(data) / factor
  data <- data[names(data) %in% seqnames]
  export.bw(data, paste0("/DATA/usr/m.trauernicht/projects/ATAC_TF_reporter_comparison/ATAC_seq_analyses/bigwig/mouse_selected_mt20260420/", group, ".bw"))
}
