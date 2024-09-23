# Library statements ------------------------------------------------------

library(here)
library(data.table)
library(dplyr, include.only = "case_when")
library(readr, include.only = "write_tsv")
library(stringr, include.only = "str_extract")
library(glue)

# Declare samples ---------------------------------------------------------

df <- tibble::tribble(
  ~basename, ~celltype, ~treatment,
  # Run 7822 + 7823 + 7824
  "U2OS_PMA", "U2OS", "PMA",
  "U2OS_Heat", "U2OS", "Heat",
  "U2OS_DMSO", "U2OS", "DMSO",
  "MCF7_Hex", "MCF7", "Hexestrol",
  "MCF7_DMSO", "MCF7", "DMSO",
  "U2OS_Calcitriol", "U2OS", "Calcitriol"
)

# Files -------------------------------------------------------------------

# List all files from the runs indicated above
runs <- paste0(sort(unique(df$basename)), collapse = "|")
gcffiles <- list.files("/DATA/usr/m.trauernicht/projects/ATAC_TF_reporter_comparison/ATAC_seq_analyses/data/human/comb", pattern = runs, recursive = TRUE, full.names = TRUE)
base_gcf <- basename(gcffiles)

# Find the matching fastq files
read1 <- gcffiles[pmatch(with(df, glue("{basename}_R1")), base_gcf)]
read1 <- read1[!is.na(read1)]
read2 <- gcffiles[pmatch(with(df, glue("{basename}_R2")), base_gcf)]
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

dir <- "/DATA/usr/m.trauernicht/projects/ATAC_TF_reporter_comparison/ATAC_seq_analyses/data/human/comb"
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
  df, "/DATA/usr/m.trauernicht/projects/ATAC_TF_reporter_comparison/ATAC_seq_analyses/ATAC_seq_samples_human_comb.tsv"
  
)
