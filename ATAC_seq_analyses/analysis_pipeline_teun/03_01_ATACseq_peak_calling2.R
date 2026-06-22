# Libraries ---------------------------------------------------------------

library(here)
library(glue)

path_bgzip <- "/DATA/usr/t.filipovska/software/Miniconda3/pkgs/tabix-0.2.6-ha92aebf_0/bin"

# Directories and data ----------------------------------------------------

files_df <- readRDS("/home/m.trauernicht/projects/ATAC_TF_reporter_comparison/ATAC_seq_analyses/rds/bamfile_atacseq_metadata_mES_selected_mt20260420.rds")
macs2 <- "/DATA/usr/t.filipovska/software/Miniconda3/envs/tf_activity/bin/macs2"
peak_dir <- "/DATA/usr/m.trauernicht/projects/ATAC_TF_reporter_comparison/ATAC_seq_analyses/bed_peaks/mouse_selected_mt20260420"
bedtools <- "/usr/bin/bedtools"
genome_sizes <- "/DATA/usr/m.trauernicht/data/genomes/mm10/mm10.chrom.sizes"

exps <- split(files_df$tabix_file, files_df$run)
#exps[names(exps) != "technical"]

# Functions ---------------------------------------------------------------

merge_tabixes <- function(files, out_file = tempfile(fileext = ".bed.gz")) {
  files <- paste0(files, collapse = " ")
  cmd <- glue("cat {files} > {out_file}")
  system(cmd)
  cmd <- glue("{path_bgzip}/bgzip -d {out_file}")
  system(cmd)
  gsub(".gz$", "", out_file)
}

call_macs2 <- function(file, name) {
  cmd <- glue(
    "{macs2} callpeak -t {file} -f BEDPE -g mm -n {name}", # change mm to hg for human
    " --nomodel --outdir {peak_dir} --keep-dup all"
  )
  system(cmd)
}

make_fixed_peaks <- function(narrowpeak_file, fixed_file) {
  cmd <- paste(
    "awk 'BEGIN{OFS=\"\\t\"} {summit=$2+$10; if(summit<0) summit=0; print $1, summit, summit+1, $4, $9, \".\"}'",
    shQuote(narrowpeak_file),
    "|",
    bedtools,
    "slop -i - -g",
    shQuote(genome_sizes),
    "-l 249 -r 250",
    "| awk 'BEGIN{OFS=\"\\t\"} {print $1,$2,$3,$4,$5}' >",
    shQuote(fixed_file)
  )
  status <- system(cmd)
  stopifnot(status == 0)
}


for (expname in names(exps)) {
  tabixes <- exps[[expname]]
  peakfile <- file.path(peak_dir, glue("{expname}_peaks.narrowPeak"))
  fixedfile <- file.path(peak_dir, glue("{expname}_peaks_fixed.bed"))
  if (!file.exists(peakfile)) {
    temp <- merge_tabixes(tabixes)
    call_macs2(temp, expname)
    unlink(temp)
  }
  if (!file.exists(fixedfile)) {
    make_fixed_peaks(peakfile, fixedfile)
  }
}

