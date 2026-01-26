# snakemake_usearch
Simple snakemake workflow based on the standard [usearch11 SOP](https://drive5.com/usearch/manual/uparse_pipeline.html).

## Brief overview and description of steps (rules)
The following steps are performed in order, see the individual files for additional details:

| File | Rule | Description |
| --- | --- | --- |
| `01-sample_prep.smk` | `sample_prep` | For each sample, decompress, concatenate all files into one, count reads before filtering, relabel reads using `usearch -fastx_relabel`, run `filtlong`, and count reads after filtering. |
| | `concatenate_total_reads_files` | Concatenates all files with number of reads before and after filtering in each sample into a single file. |
| `02-denoise.smk` | `concat_all` | Concatenates all filtered reads from all samples into a single file to be able to generate ASVs/zOTUs. |
| | `trim_primers` | Trim primers and orient reads using `cutadapt`. |
| | `derep` | Dereplicate reads using `usearch -fastx_uniques`. |
| | `unoise` | Generates ASVs/zOTUs using `usearch -unoise3`. |
| `04-sintax.smk` | `sintax` | Predict taxonomy of the ASVs/zOTUs using `usearch -sintax`. |
| `05-abund_table.smk` | `abund_table` | For each sample, estimate ASV/zOTU abundances in each sample by mapping the raw, unfiltered reads against the ASVs/zOTUs using `usearch -otutab`. This is MUCH faster to do in parallel and merge afterwards compared to running a single `usearch -otutab` command, which doesn't scale linearly with more threads. |
| | `merge_abund_tables` | Merge all abundance tables into a single table using `usearch -otutab_merge`. |
| | `rarefy_abund_table` | (optional) Rarefy abundance table using `usearch -otutab_rare`. |

## Usage
First install snakemake and the required software into a conda environment or use the container as described below. Then run using fx:
```
conda activate snakemake_usearch
snakemake --jobs 96
```

Use an executor if you are running on a HPC cluster. See the `slurm_submit.sbatch` for an example when running on a SLURM cluster.

## Requirements
Install the required software by using either the provided `Dockerfile` or `environment.yml` file to build a Docker container or conda environment with all the required tools, see below.

### Docker
The pre-built Docker container available from [`ghcr.io/kasperskytte/snakemake_usearch:main`](https://github.com/KasperSkytte/snakemake_usearch/pkgs/container/snakemake_usearch) is built from the `Dockerfile` and includes all required software.

### Conda
Requirements are listed in `environment.yml`. To create as a conda environment simply run:
```
conda env create --file environment.yml -n snakemake_usearch
```
