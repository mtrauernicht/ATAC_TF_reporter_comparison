# Systematic comparison of estimates of transcription factor activity by ATAC-seq and multiplexed reporter assays

## Overview

This repository contains the analysis and code for comparing transcription factor (TF) activity measurements obtained via two methods: **TF Reporter Assays** and indirect detection using **ATAC-seq**.

## Goals of the Project

- To evaluate the correlation between TF activity measured by **reporter assays** and **ATAC-seq**.
- To identify cases where one method provides better sensitivity or specificity than the other.
- To provide insights into how well TF reporter assays reflect the actual transcriptional and chromatin state changes that occur in response to TF perturbation.

## Project Structure

- **`/data/`**: Contains datasets used for the analysis, including:
  - TF reporter assay data
  - ATAC-seq data
  - RNA-seq data
- **`/ATAC_seq_analyses/`**: Contains all analysis pipelines to analyze the raw ATAC-seq data.
- **`/TF_reporter_analyses/`**: Contains computation of TF activity from TF reporter assays + computation of TF activity from ATAC-seq & RNA-seq and their comparison.


## Contact

For any questions, feel free to reach out to [m.trauernicht@nki.nl].

---

