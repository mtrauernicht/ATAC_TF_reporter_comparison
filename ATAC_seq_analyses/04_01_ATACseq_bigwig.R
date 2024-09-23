# Libraries ---------------------------------------------------------------

library(here)
library(glue)
library(data.table)

library(GenomicRanges)
library(rtracklayer)
library(IRanges)
library(GenomeInfoDb)


# Directories and data ----------------------------------------------------

files_df <- readRDS(here("ATAC_seq_analyses/rds", "bamfile_atacseq_metadata.rds"))
bw_dir <- "/DATA/usr/m.trauernicht/projects/ATAC_TF_reporter_comparison/ATAC_seq_analyses/bigwig/7784"

# Split up pools ----------------------------------------------------------

files_df <- files_df[files_df$run == "7784",]

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
  
  n <- sum(width(data))
  factor <- gsize / n
  
  data <- coverage(data) * factor
  data <- data[names(data) %in% seqnames]
  export.bw(data, paste0("/DATA/usr/m.trauernicht/projects/ATAC_TF_reporter_comparison/ATAC_seq_analyses/bigwig/", group, ".bw"))
}
