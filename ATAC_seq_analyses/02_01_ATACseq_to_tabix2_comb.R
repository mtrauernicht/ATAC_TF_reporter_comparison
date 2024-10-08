# Libraries ---------------------------------------------------------------
library(GenomicAlignments)
library(GenomicRanges)
library(Rsamtools)
library(S4Vectors)
library(data.table, exclude = c("first", "second", "shift", "last"))
library(here)
library(glue)

path_bgzip <- "/DATA/usr/t.filipovska/software/Miniconda3/pkgs/tabix-0.2.6-ha92aebf_0/bin"

# Files and directories ---------------------------------------------------
files <- readRDS(here("ATAC_seq_analyses/rds", "bamfile_atacseq_metadata.rds"))
dir <- "/DATA/usr/m.trauernicht/projects/ATAC_TF_reporter_comparison/ATAC_seq_analyses/"
dir_in  <- glue("{dir}bam/7784")
dir_out <- glue("{dir}tabix/7784_comb")

# Merge replicate BAM files -----------------------------------------------
# Assuming 'files' has a column 'treatment' to identify replicates
files$treatment <- as.factor(files$treatment)  # Ensure 'treatment' is a factor

# Create a list of BAM files grouped by treatment
bamfiles_list <- split(files$bam_file, files$treatment)

# Function to merge BAM files
merge_bam_files <- function(bamfiles, output_file) {
  if (file.exists(output_file)) {
    file.remove(output_file)
  }
  mergeBam(bamfiles, destination = output_file, indexDestination = TRUE)
}

# Merged BAM files
merged_bam_files <- sapply(names(bamfiles_list), function(treatment) {
  bamfiles <- bamfiles_list[[treatment]]
  output_file <- glue("{dir_in}/{treatment}_merged.bam")
  merge_bam_files(bamfiles, output_file)
  return(output_file)
})

# Create BamFileList
bamfiles <- BamFileList(merged_bam_files, yieldSize = 1e6L, asMates = TRUE)

# Blacklist from https://dozmorovlab.github.io/HiCcompare/reference/hg38_blacklist.html
#blacklist <- fread("/DATA/usr/t.filipovska/projects/ATAC-seq_TF-footprinting/analysis/references/hg38.blacklist.bed.gz")
blacklist <- fread("/DATA/usr/t.filipovska/projects/ATAC-seq_TF-footprinting/analysis/references/ENCFF547MET.bed.gz")
blacklist <- with(blacklist, GRanges(V1, IRanges(V2, V3)))
blacklist <- as(blacklist, "GNCList")

# Setup parameters --------------------------------------------------------
setDTthreads(1)

param <- ScanBamParam(
  flag = scanBamFlag(
    isPaired = TRUE, isProperPair = TRUE, isSecondaryAlignment = FALSE,
    isUnmappedQuery = FALSE, hasUnmappedMate = FALSE
  ),
  mapqFilter = 10L
)

# Functions ---------------------------------------------------------------
shift_fragments <- function(frags, positive = 4, negative = -5) {
  read1 <- tn5shift(GenomicAlignments::first(frags), positive, negative)
  read2 <- tn5shift(GenomicAlignments::second(frags), positive, negative)
  punion(read1, read2, fill.gap = TRUE, ignore.strand = TRUE)
}

tn5shift <- function(x, positive = 4, negative = -5) {
  x <- resize(granges(x), fix = "start", 1)
  offset <- strand(x) == "+"
  offset <- Rle(ifelse(runValue(offset), positive, negative), runLength(offset))
  GenomicRanges::shift(x, decode(offset))
}

write_bed <- function(granges, file) {
  granges <- data.table(
    chrom = decode(seqnames(granges)),
    start = start(granges),
    end = end(granges)
  )
  fwrite(granges, file, append = TRUE, quote = FALSE, sep = "\t", col.names = FALSE)
  file
}

finalise_tabix <- function(file) {
  sorter <- glue("sort -k1,1 -k2,2n -o {file} {file}")
  system(sorter)
  zipper <- glue("{path_bgzip}/bgzip {file}")
  system(zipper)
  file <- paste0(file, ".gz")
  index <- indexTabix(file, seq = 1, start = 2, end = 3)
  TabixFile(file, index)
}

read_bam <- function(file) {
  readGAlignmentPairs(file, param = param)
}

# Convert -----------------------------------------------------------------
bedfiles <- path(bamfiles)
bedfiles <- paste0(dir_out, "/", basename(bedfiles))
bedfiles <- gsub(".bam$", ".bed", bedfiles)
tabixfiles <- paste0(bedfiles, ".gz")

for (i in seq_along(bamfiles)) {
  bamfile <- bamfiles[[i]]
  outfile <- bedfiles[[i]]

  if (file.exists(tabixfiles[[i]])) {
    next
  }

  open(bamfile)
  while (length({data <- read_bam(bamfile)})) {
    data <- shift_fragments(data)
    data <- data[width(data) > 19]
    data <- data[!overlapsAny(data, blacklist)]
    outfile <- write_bed(data, outfile)
  }
  close(bamfile)

  tbx <- finalise_tabix(outfile)
}

files$tabix_file <- tabixfiles
saveRDS(files, here("ATAC_seq_analyses/rds", "bamfile_atacseq_metadata2_comb.rds"))

