# snakemake_usearch
Simple snakemake workflow based on the standard [usearch11 SOP](https://drive5.com/usearch/manual/uparse_pipeline.html).


## Brief overview and description of steps (rules)

The Snakemake workflow is modular and organized into sequential preprocessing steps followed by ASV inference and downstream analyses. Individual rules are defined in separate `.smk` files; see the corresponding rule files for implementation details.

### Global preprocessing steps

| File | Rule(s) | Description |
| --- | --- | --- |
| `rules/01-sample_prep.smk` | `sample_prep` | Decompress and merge demultiplexed FASTQ files per sample as an initial preprocessing step prior to downstream analyses. |
| `rules/02-primertrim.smk` | `primertrim` | Trim primer sequences and orient reads using **cutadapt**. |
| `rules/03-subsample.smk` | `subsample` | Subsample reads into user-defined subsets specified in the configuration file (e.g. to evaluate sequencing depth effects). |

---

### Per-subset ASV processing

The following steps are executed **independently for each subsampled read set**:

| File | Rule(s) | Description |
| --- | --- | --- |
| `rules/04-denoise.smk` | `denoise` | Perform ASV inference using **USEARCH `-unoise3`** on each subset. |
| `rules/05-sintax.smk` | `sintax` | Assign taxonomy to ASVs using **USEARCH `-sintax`** with a user-specified reference database. |
| `rules/06-append_asv_counts.smk` | `append_asv_counts` | (Optional) Match inferred ASVs against a user-defined reference ASV set and append count information; this step can be skipped if not required. |
| `rules/07-abund_table.smk` | `abund_table` | Estimate ASV abundances per sample by mapping reads to ASVs using **USEARCH `-otutab`**. |

---


## Usage
First install snakemake and the required software into a conda environment:

```
conda activate snakemake_usearch
snakemake --jobs 96
```

Use an executor if you are running on a HPC cluster. See the `slurm_submit.sbatch` for an example when running on a SLURM cluster.

## Requirements
Install the required software by using the `environment.yml` file to build a conda environment with all the required tools, see below.

Requirements are listed in `environment.yml`. To create as a conda environment simply run:

```
conda env create --file environment.yml -n snakemake_usearch
```
