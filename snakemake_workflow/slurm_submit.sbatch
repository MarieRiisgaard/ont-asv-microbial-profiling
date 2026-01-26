#!/usr/bin/bash -l
#SBATCH --job-name=snakemake_usearch
#SBATCH --output=job_%j_%x.out
#SBATCH --cpus-per-task=1
#SBATCH --time=0-06:00:00
#SBATCH --mem=1G

# Exit on first error and if any variables are unset
set -eu

# Activate the conda environment
conda activate snakemake_usearch

# Optionally clear files from a previous run before starting
#rm -rf logs tmp output
#snakemake --unlock

# Run snakemake
snakemake --executor slurm --jobs 96
