# Libraries ---------------------------------------------------------------

library(SummarizedExperiment)
library(GenomicRanges)
library(here)
library(rtracklayer)
library(data.table)

# Data --------------------------------------------------------------------

meta <- readRDS("/DATA/usr/m.trauernicht/projects/ATAC_TF_reporter_comparison/ATAC_seq_analyses/rds/bamfile_atacseq_metadata2_comb.rds")

peaks_7784 <- import("/DATA/usr/m.trauernicht/projects/ATAC_TF_reporter_comparison/ATAC_seq_analyses/bed_peaks/20240626_comb/7784_peaks.narrowPeak")
peaks_7822 <- import("/DATA/usr/m.trauernicht/projects/ATAC_TF_reporter_comparison/ATAC_seq_analyses/bed_peaks/20240626_comb/7822_peaks.narrowPeak")
peaks <- c(peaks_7784, peaks_7822)
peaks <- as(peaks, "GNCList")


# Metadata formatting -----------------------------------------------------

#meta <- meta[meta$run == "7784", ]
meta <- meta[order(interaction(meta$celltype, meta$treatment)),]
#rep <- tapply(meta$celltype, interaction(meta$celltype, meta$treatment, meta$replicate))
#meta$rep <- unlist(rep)
meta$alias <- paste(meta$celltype, meta$treatment, sep = "_")
meta$celltype <- factor(meta$celltype, levels = c("mES"))
meta$treatment <- factor(meta$treatment, levels = c("2i_LIF", "FOXA1_OE", "FOXA1_ctrl", "SOX2_DEG", "SOX2_ctrl", "POU5F1_DEG", "POU5F1_ctrl",
                                                    "LIF", "CH", "PD"))
#meta$replicate <- factor(meta$replicate, levels=c("R1", "R2", "R3"))

# Count -------------------------------------------------------------------

mat <- vapply(meta$tabix_file, function(file) {
  dat <- data.table::fread(file)
  dat <- with(dat, GRanges(V1, IRanges(V2, V3)))
  olaps <- findOverlaps(dat, peaks)
  olaps <- data.table(olap = to(olaps))
  olaps <- olaps[, .N, by = olap]
  olaps <- olaps[.(seq_along(peaks)), on = "olap"]
  return(olaps$N)
}, integer(NROW(peaks)))
mat[is.na(mat)] <- 0
colnames(mat) <- meta$alias


# Summarise Experiment ----------------------------------------------------

exp <- SummarizedExperiment(
  assays = mat,
  rowRanges = as(peaks, "GRanges"),
  colData = meta
)

saveRDS(exp, "/DATA/usr/m.trauernicht/projects/ATAC_TF_reporter_comparison/ATAC_seq_analyses/rds/experiment_7784_20240724_comb.rds")
