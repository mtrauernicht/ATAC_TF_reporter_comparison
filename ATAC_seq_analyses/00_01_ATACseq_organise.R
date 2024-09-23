# Library statements ------------------------------------------------------

library(here)
library(data.table)
library(dplyr, include.only = "case_when")
library(readr, include.only = "write_tsv")
library(stringr, include.only = "str_extract")
library(glue)

# Declare samples ---------------------------------------------------------

df <- tibble::tribble(
  ~run, ~basename, ~celltype, ~treatment, ~replicate,
  # Run 7784 + 7822
  "7784", "10_mES_FOXA1-OE_ctrl_r1_CGTACTAG-AGAGGATA_S10", "mES", "FOXA1_ctrl", "R2",
  "7784", "11_mES_2i_r2_TAAGGCGA-CTCCTTAC_S11", "mES", "LIF", "R3",
  "7784", "12_mES_LIF_PD_r2_CGTACTAG-CTCCTTAC_S12", "mES", "CH", "R3",
  "7784", "13_mES_LIF_CH_r2_AGGCAGAA-CTCCTTAC_S13", "mES", "PD", "R3", 
  "7784", "14_mES_2i_LIF_r2_TCCTGAGC-CTCCTTAC_S14", "mES","2i_LIF", "R3", 
  "7784", "15_mES_POU5F1_DEG_r2_GGACTCCT-CTCCTTAC_S15", "mES", "POU5F1_DEG", "R3", 
  "7784", "16_mES_POU5F1_ctrl_r2_TAGGCATG-CTCCTTAC_S16", "mES", "POU5F1_ctrl", "R3",
  "7784", "17_mES_SOX2_DEG_r2_CTCTCTAC-CTCCTTAC_S17", "mES", "SOX2_DEG", "R3",
  "7784", "18_mES_SOX2_ctrl_r2_CAGAGAGG-CTCCTTAC_S18", "mES", "SOX2_ctrl", "R3",
  "7784", "19_mES_FOXA1-OE_r2_TAAGGCGA-TATGCAGT_S19", "mES","FOXA1_OE", "R3",
  "7784", "1_mES_2i_r1_TAAGGCGA-ATAGAGAG_S1", "mES","LIF", "R2",
  "7784", "20_mES_FOXA1-OE_ctrl_r2_CGTACTAG-TATGCAGT_S20", "mES", "FOXA1_ctrl", "R3",
  "7784", "2_mES_LIF_PD_r1_CGTACTAG-ATAGAGAG_S2", "mES", "CH", "R2",
  "7784", "3_mES_LIF_CH_r1_AGGCAGAA-ATAGAGAG_S3", "mES", "PD", "R2",
  "7784", "4_mES_2i_LIF_r1_TCCTGAGC-ATAGAGAG_S4", "mES", "2i_LIF", "R2",
  "7784", "5_mES_POU5F1_DEG_r1_GGACTCCT-ATAGAGAG_S5", "mES","POU5F1_DEG", "R2",
  "7784", "6_mES_POU5F1_ctrl_r1_TAGGCATG-ATAGAGAG_S6", "mES","POU5F1_ctrl", "R2",
  "7784", "7_mES_SOX2_DEG_r1_CTCTCTAC-ATAGAGAG_S7", "mES", "SOX2_DEG", "R2",
  "7784", "8_mES_SOX2_ctrl_r1_CAGAGAGG-ATAGAGAG_S8", "mES","SOX2_ctrl", "R2",
  "7784", "9_mES_FOXA1-OE_r1_TAAGGCGA-AGAGGATA_S9", "mES", "FOXA1_OE", "R2",
  "7822", "10_mES_SOX2_ctrl_r1_TAAGGCGA-AGAGGATA_S10", "mES", "SOX2_ctrl", "R1",
  "7822", "11_mES_FOXA1_OE_r1_CGTACTAG-AGAGGATA_S11", "mES", "FOXA1_OE", "R1",
  "7822", "12_mES_FOXA1_ctrl_r1_AGGCAGAA-AGAGGATA_S12", "mES", "FOXA1_ctrl", "R1",
  "7822", "3_mES_2i_r1_TAAGGCGA-GCGATCTA_S3", "mES", "LIF", "R1",
  "7822", "4_mES_LIF_PD_r1_AGGCAGAA-ATAGAGAG_S4", "mES", "CH", "R1",
  "7822", "5_mES_LIF_CH_r1_TCCTGAGC-ATAGAGAG_S5", "mES", "PD", "R1",
  "7822", "6_mES_2i_LIF_r1_GGACTCCT-ATAGAGAG_S6", "mES", "2i_LIF", "R1",
  "7822", "7_mES_POU5F1_DEG_r1_TAGGCATG-ATAGAGAG_S7", "mES", "POU5F1_DEG", "R1",
  "7822", "8_mES_POU5F1_ctrl_r1_CTCTCTAC-ATAGAGAG_S8", "mES", "POU5F1_ctrl", "R1",
  "7822", "9_mES_SOX2_DEG_r1_CAGAGAGG-ATAGAGAG_S9", "mES", "SOX2_DEG", "R1"
)

# Files -------------------------------------------------------------------

# List all files from the runs indicated above
runs <- paste0(sort(unique(df$run)), collapse = "|")
gcffiles <- list.files("/shared/gcf", pattern = runs, recursive = TRUE, full.names = TRUE)
base_gcf <- basename(gcffiles)

# Find the matching fastq files
read1 <- gcffiles[pmatch(with(df, glue("{run}_{basename}_R1")), base_gcf)]
read1 <- read1[!is.na(read1)]
read2 <- gcffiles[pmatch(with(df, glue("{run}_{basename}_R2")), base_gcf)]
read2 <- read2[!is.na(read2)]

# Check if all samples have been found
check <- any(is.na(read1) | is.na(read2))
stopifnot(!check)

# Check filename validity -------------------------------------------------

names <- c(read1, read2)
names <- gsub("R1|R2", "", names)
names <- matrix(names, ncol = 2)

# If removing R1 and R2, all parallel filenames should be equal
check <- all(names[, 1] == names[, 2])
stopifnot(check)

# Symlinking files --------------------------------------------------------

dir <- "/DATA/usr/m.trauernicht/projects/ATAC_TF_reporter_comparison/ATAC_seq_analyses/data"
new_read1 <- glue("{dir}/{basename(read1)}")
new_read2 <- glue("{dir}/{basename(read2)}")

file_is_linked <- function(x) {
  f <- list.files(dirname(x)[1], full.names = TRUE)
  x %in% f
}

check <- file_is_linked(new_read1) & file_is_linked(new_read2)

cmd_r1 <- glue("ln -s {read1} {new_read1}")
cmd_r2 <- glue("ln -s {read2} {new_read2}")

for (i in seq_along(check)[!check]) {
  system(cmd_r1[i])
  system(cmd_r2[i])
}

# Update data with fastq locations ----------------------------------------

df <- transform(
  df,
  read1 = new_read1,
  read2 = new_read2
)

# Export ------------------------------------------------------------------

write_tsv(
  df, "/DATA/usr/m.trauernicht/projects/ATAC_TF_reporter_comparison/ATAC_seq_analyses/ATAC_seq_samples.tsv"
  
)
