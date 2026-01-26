# ONT ASV-resolved Microbial Profiling

This repository contains the scripts and Snakemake workflows used in the study:

**_Improved Nanopore sequencing accuracy enables ASV-resolved microbial profiling_**

The repository implements a reproducible workflow for amplicon sequence variant (ASV) inference from Oxford Nanopore Technologies (ONT) 16S rRNA gene amplicon data, including preprocessing, ASV generation, and downstream analyses.

---

## Overview of the workflow

**ONT amplicon preprocessing**

The snakemake pipeline is also availible: https://github.com/MarieRiisgaard/snakemake_usearch_primertrim_and_subset where it can be deployed directly

Input are demultiplexed reads and work for both ONT and PacBio data. 

The Snakemake workflow performs:
   - Primer trimming using CutAdapt using primers specified by the user
   - Read subsetting (based on user specified input)
   - High-accuracy ASV generation from ONT amplicons
   - Length filtering
   - Classification of the generated ASVs 

## Downstream analyses

### Scripts

The following helper scripts are used for specific preprocessing, ASV generation, and downstream steps:

- **`calculate_mean_qscore.sh`**  
  Calculate mean Phred Q-scores for the subsampled reads, used to summarize read-quality across sample types and sequencing platforms.

- **`chimera_check_zymo.sh`**  
  Screen generated amplicons/ASVs for chimeras using *USEARCH `uchime`* against a reference database (ZymoMock reference).

- **`generate_in_silico_amplicons.sh`**  
  Generate in silico reference amplicons by trimming/cutting reference sequences using the relevant primer sets (used for constructing ZymoMock in silico references).

- **`make_global_asvs_search_ref.sh`**  
  Merge all generated ASVs per primer set across all subsets and samples, for both ONT and PacBio, to create a global ASV reference set. This global set provides stable IDs for identifying identical ASVs across datasets.

- **`Parallel_ASVtab_OTUtab_mapping.sh`**  
  Construct an OTU/ASV table by mapping reads to a user-specified ASV set (e.g., the global ASV reference).

- **`process_DADA2_amplicons.sh`**  
  Postprocess amplicons generated with DADA2 and preparing outputs for downstream comparisons.

- **`qiime2_ASVs_generation.sh`**  
  Run the *QIIME 2* workflow for ASV generation and core QIIME 2 processing steps.

- **`qiime2_slurm_submit.sh`**  
  Submit the QIIME 2 workflow to an HPC/SLURM environment (job wrapper for cluster execution).

- **`article_script_final.Rmd`**  
  Rmarkdown to generate figures in the article.

---

## Requirements

The pipeline is designed to run on Linux systems and HPC environments.

### Core dependencies
- Snakemake (≥ 7)
- Conda / Mamba
- Python (via Conda)
- R (for postprocessing)

All software dependencies for the main workflow are defined in: snakemake_workflow/environment.yml

---


## Data availability

Raw sequencing data are available via the Sequence Read Archive (SRA) under the accession number(s) reported in the manuscript.

---

## Citation

If you use this workflow, please cite:

Riisgaard-Jensen, M. et al.
Improved Nanopore sequencing accuracy enables ASV-resolved microbial profiling
[Journal, Year]


