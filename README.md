# Systematic comparison of estimates of transcription factor activity by ATAC-seq and multiplexed reporter assays

<a href="https://doi.org/10.5281/zenodo.18019661"><img src="https://zenodo.org/badge/819850868.svg" alt="DOI"></a>

## Overview

This repository contains the analysis and code for comparing transcription factor (TF) activity measurements obtained via two methods: **TF Reporter Assays** and indirect detection using **ATAC-seq**.

## Goals of the Project

- To evaluate the correlation between TF activity measured by **reporter assays** and **ATAC-seq**.
- To identify cases where one method provides better sensitivity or specificity than the other.
- To provide insights into how well TF reporter assays reflect the actual transcriptional and chromatin state changes that occur in response to TF perturbation.

## Project Structure

- **`ATAC_reporter_TF_activity_comparison.Rmd`**: Contains all code used to make Figures. 
- **`/ATAC_seq_analyses/`**: Contains all analysis pipelines to analyze the raw ATAC-seq data and the analysis using chromVAR.
- **`/TF_reporter_analyses/`**: Contains computation of TF activity from TF reporter assays using primetime.


## Contact

For any questions, feel free to reach out to [m.trauernicht@nki.nl].

---

