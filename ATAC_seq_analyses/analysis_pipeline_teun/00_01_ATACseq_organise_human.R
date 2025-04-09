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
  # Run 7822 + 7823 + 7824
  "7822", "1_U2OS_PMA_r1_TAAGGCGA-ATAGAGAG_S1", "U2OS", "PMA", "R1",
  "7822", "2_U2OS_Heat_r1_CGTACTAG-ATAGAGAG_S2", "U2OS", "Heat", "R1",
  "7823", "10_U2OS_DMSO_r3_TAGGCATG-TCTACTCT_S10", "U2OS", "DMSO", "R3",
  "7823", "11_MCF7_Hex_r3_CTCTCTAC-TCTACTCT_S11", "MCF7", "Hexestrol", "R3",
  "7823", "12_MCF7_DMSO_r3_CAGAGAGG-TCTACTCT_S12", "MCF7", "DMSO", "R3",
  "7823", "1_U2OS_Calcitriol_r1_TAAGGCGA-GCGATCTA_S1", "U2OS", "Calcitriol", "R1",
  "7823", "2_U2OS_Calcitriol_r2_AGGCAGAA-GCGATCTA_S2", "U2OS", "Calcitriol", "R2",
  "7823", "3_U2OS_PMA_r2_TCCTGAGC-GCGATCTA_S3", "U2OS", "PMA", "R2",
  "7823", "4_U2OS_Heat_r2_GGACTCCT-GCGATCTA_S4", "U2OS", "Heat", "R2",
  "7823", "5_U2OS_DMSO_r2_TAGGCATG-GCGATCTA_S5", "U2OS", "DMSO", "R2",
  "7823", "6_MCF7_Hex_r2_CTCTCTAC-GCGATCTA_S6", "MCF7", "Hexestrol", "R2",
  "7823", "7_MCF7_DMSO_r2_CAGAGAGG-GCGATCTA_S7", "MCF7", "DMSO", "R2",
  "7823", "8_U2OS_Calcitriol_r3_AGGCAGAA-TCTACTCT_S8", "U2OS", "Calcitriol", "R3",
  "7823", "9_U2OS_PMA_r3_TCCTGAGC-TCTACTCT_S9", "U2OS", "PMA", "R3",
  "7824", "17_ATAC_U2OS_DMSO_r1_TAAGGCGA-GCGATCTA_S17", "U2OS", "DMSO", "R1"
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

dir <- "/DATA/usr/m.trauernicht/projects/ATAC_TF_reporter_comparison/ATAC_seq_analyses/data/human"
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
  df, "/DATA/usr/m.trauernicht/projects/ATAC_TF_reporter_comparison/ATAC_seq_analyses/ATAC_seq_samples_human.tsv"
  
)
