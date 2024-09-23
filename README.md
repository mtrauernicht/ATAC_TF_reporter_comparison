# Comparison of multiplexed TF activity detection methods

## Overview

This repository contains the analysis and code for comparing transcription factor (TF) activity measurements obtained via two methods: **TF Reporter Assays** and indirect detection using **ATAC-seq** and **RNA-seq**. The goal is to systematically evaluate how well TF reporter assays correlate with the more indirect approaches of chromatin accessibility (ATAC-seq) and transcriptomic changes (RNA-seq) to measure TF activity.

## Goals of the Project

- To evaluate the correlation between **direct** TF activity measurement (reporter assays) and **indirect** methods (ATAC-seq and RNA-seq).
- To identify cases where one method provides better sensitivity or specificity than the other.
- To provide insights into how well TF reporter assays reflect the actual transcriptional and chromatin state changes that occur in response to TF perturbation.

## Project Structure

- **`/data/`**: Contains datasets used for the analysis, including:
  - TF reporter assay data
  - ATAC-seq data
  - RNA-seq data
- **`/ATAC_seq_analyses/`**: Contains all analysis pipelines to analyze the raw ATAC-seq data.
- **`/RNA_seq_analyses/`**: Contains analyses to process the RNA-seq data.
- **`/TF_reporter_analyses/`**: Contains computation of TF activity from TF reporter assays + computation of TF activity from ATAC-seq & RNA-seq and their comparison.
- **`/docs/`**: Documentation files explaining the methodology, findings, and supplementary information.


## Contact

For any questions, feel free to reach out to [m.trauernicht@nki.nl].

---

