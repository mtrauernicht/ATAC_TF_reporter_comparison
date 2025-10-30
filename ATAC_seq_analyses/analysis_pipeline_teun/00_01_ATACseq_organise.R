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
"7822",	"6_mES_2i_LIF_r1_GGACTCCT-ATAGAGAG_S6",	"mES", "2i_LIF",	"R1",
"7822",	"4_mES_LIF_PD_r1_AGGCAGAA-ATAGAGAG_S4",	"mES", "CH",	"R1",
"7822",	"3_mES_2i_r1_TAAGGCGA-GCGATCTA_S3",	"mES", "LIF",	"R1",
"7822",	"8_mES_POU5F1_ctrl_r1_CTCTCTAC-ATAGAGAG_S8",	"mES", "POU5F1_ctrl",	"R1",
"7822",	"7_mES_POU5F1_DEG_r1_TAGGCATG-ATAGAGAG_S7",	"mES", "POU5F1_DEG",	"R1",
"7822",	"10_mES_SOX2_ctrl_r1_TAAGGCGA-AGAGGATA_S10",	"mES", "SOX2_ctrl",	"R1",
"7822",	"9_mES_SOX2_DEG_r1_CAGAGAGG-ATAGAGAG_S9",	"mES", "SOX2_DEG",	"R1",
"7784",	"4_mES_2i_LIF_r1_TCCTGAGC-ATAGAGAG_S4",	"mES", "2i_LIF",	"R2",
"7784",	"2_mES_LIF_PD_r1_CGTACTAG-ATAGAGAG_S2",	"mES", "CH",	"R2",
"8010",	"1_mES_FBS_r2_TAAGGCGA-TACTCCTT_S1",	"mES", "FBS",	"R2",
"8010",	"2_mES_FK_r2_CGTACTAG-TACTCCTT_S2",	"mES", "FK",	"R2",
"8010",	"18_mES_Heat_r2_CTCTCTAC-TATGCAGT_S18",	"mES", "Heat",	"R2",
"7784",	"1_mES_2i_r1_TAAGGCGA-ATAGAGAG_S1",	"mES", "LIF",	"R2",
"8010",	"12_mES_LIF_1_r2_TAAGGCGA-TATGCAGT_S12",	"mES", "LIF_1h",	"R2",
"8010",	"11_mES_LIF_3_r2_CTCTCTAC-CTCCTTAC_S11",	"mES", "LIF_3h",	"R2",
"8010",	"10_mES_LIF_6_r2_TAGGCATG-CTCCTTAC_S10",	"mES", "LIF_6h",	"R2",
"8010",	"9_mES_NT_r2_GGACTCCT-CTCCTTAC_S9",	"mES", "NT",	"R2",
"7784",	"6_mES_POU5F1_ctrl_r1_TAGGCATG-ATAGAGAG_S6",	"mES", "POU5F1_ctrl",	"R2",
"7784",	"5_mES_POU5F1_DEG_r1_GGACTCCT-ATAGAGAG_S5",	"mES", "POU5F1_DEG",	"R2",
"7784",	"8_mES_SOX2_ctrl_r1_CAGAGAGG-ATAGAGAG_S8",	"mES", "SOX2_ctrl",	"R2",
"7784",	"7_mES_SOX2_DEG_r1_CTCTCTAC-ATAGAGAG_S7",	"mES", "SOX2_DEG",	"R2",
"8010",	"8_mES_TFCP2L1_r2_TCCTGAGC-CTCCTTAC_S8",	"mES", "TFCP2L1",	"R2",
"7784",	"14_mES_2i_LIF_r2_TCCTGAGC-CTCCTTAC_S14",	"mES", "2i_LIF",	"R3",
"7784",	"12_mES_LIF_PD_r2_CGTACTAG-CTCCTTAC_S12",	"mES", "CH",	"R3",
"8011",	"17_NPC_r3_TAAGGCGA-AGAGGATA_S17",	"NPC", "DMSO",	"R3",
"8011",	"2_mES_FBS_r3_CGTACTAG-AGGCTTAG_S2",	"mES", "FBS",	"R3",
"8011",	"3_mES_FK_r3_AGGCAGAA-AGGCTTAG_S3",	"mES", "FK",	"R3",
"7784",	"20_mES_FOXA1-OE_ctrl_r2_CGTACTAG-TATGCAGT_S20",	"mES", "FOXA1_ctrl",	"R3",
"7784",	"19_mES_FOXA1-OE_r2_TAAGGCGA-TATGCAGT_S19",	"mES", "FOXA1_OE",	"R3",
"8011",	"20_mES_Heat_r3_TCCTGAGC-AGAGGATA_S20",	"mES", "Heat",	"R3",
"7784",	"11_mES_2i_r2_TAAGGCGA-CTCCTTAC_S11",	"mES", "LIF",	"R3",
"8011",	"13_mES_LIF_1_r3_GGACTCCT-ATAGAGAG_S13",	"mES", "LIF_1h",	"R3",
"8011",	"12_mES_LIF_3_r3_TCCTGAGC-ATAGAGAG_S12",	"mES", "LIF_3h",	"R3",
"8011",	"11_mES_LIF_6_r3_AGGCAGAA-ATAGAGAG_S11",	"mES", "LIF_6h",	"R3",
"8011",	"10_mES_NT_r3_CGTACTAG-ATAGAGAG_S10",	"mES", "NT",	"R3",
"7784",	"16_mES_POU5F1_ctrl_r2_TAGGCATG-CTCCTTAC_S16",	"mES", "POU5F1_ctrl",	"R3",
"7784",	"15_mES_POU5F1_DEG_r2_GGACTCCT-CTCCTTAC_S15",	"mES", "POU5F1_DEG",	"R3",
"7784",	"18_mES_SOX2_ctrl_r2_CAGAGAGG-CTCCTTAC_S18",	"mES", "SOX2_ctrl",	"R3",
"7784",	"17_mES_SOX2_DEG_r2_CTCTCTAC-CTCCTTAC_S17",	"mES", "SOX2_DEG",	"R3",
"8011",	"9_mES_TFCP2L1_r3_TAAGGCGA-ATAGAGAG_S9",	"mES", "TFCP2L1",	"R3",
"8010",	"19_NPC_r4_GGACTCCT-TACTCCTT_S19",	"NPC", "DMSO",	"R4"
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

# Copy found files to data folder --------------------------------------------
output_dir <- "/DATA/shared/projects/TF_REPORTERS/ATAC_TF_reporter_comparison/ATAC_seq_analyses/fastq_downloads"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
for (f in c(read1, read2)) {
  file.copy(f, output_dir, overwrite = FALSE)
}

# Symlinking files --------------------------------------------------------

dir <- "/DATA/shared/projects/TF_REPORTERS/ATAC_TF_reporter_comparison/ATAC_seq_analyses/data"
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
  df, "~/projects/ATAC_TF_reporter_comparison/ATAC_seq_analyses/ATAC_seq_samples.tsv"
  
)
