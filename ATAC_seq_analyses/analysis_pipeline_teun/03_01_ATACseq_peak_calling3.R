# Libraries ---------------------------------------------------------------

library(here)
library(glue)

path_bgzip <- "/DATA/usr/t.filipovska/software/Miniconda3/pkgs/tabix-0.2.6-ha92aebf_0/bin"

# Directories and data ----------------------------------------------------

files_df <- readRDS("/DATA/usr/m.trauernicht/projects/ATAC_TF_reporter_comparison/ATAC_seq_analyses/rds/bamfile_atacseq_metadata_mES.rds")
macs2 <- "/DATA/usr/t.filipovska/software/Miniconda3/envs/tf_activity/bin/macs2"
peak_dir <- "/DATA/usr/m.trauernicht/projects/ATAC_TF_reporter_comparison/ATAC_seq_analyses/bed_peaks/mouse_indiv"

# Functions ---------------------------------------------------------------

call_macs2 <- function(file, name) {
  cmd <- glue(
    "{macs2} callpeak -t {file} -f BEDPE -g mm -n {name}", # change mm to hg for human
    " --nomodel --outdir {peak_dir} --keep-dup all"
  )
  system(cmd)
}

# Iterate over each sample
for (i in seq_len(nrow(files_df))) {
  sample_file <- files_df$tabix_file[i]
  sample_name <- files_df$sample_name[i] # Assuming there's a column 'sample_name' for unique sample names
  peakfile <- paste0(peak_dir, "/", sample_name, "_peaks.narrowPeak")
  
  if (file.exists(peakfile)) {
    next
  }
  
  call_macs2(sample_file, sample_name)
}