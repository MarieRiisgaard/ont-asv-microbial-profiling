#!/bin/bash
#SBATCH --job-name=snakemake_usearch
#SBATCH --output=/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/6_snakemake_usearch_run_all/DADA2_ASV_generation/slurm_out/job_%j_%x.out
#SBATCH --cpus-per-task=10
#SBATCH --time=0-08:00:00
#SBATCH --mem=32G
#SBATCH --account=phn


cd /home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon

#############################################
# USER INPUTS — CHANGE ONLY THESE
#############################################


export basecalling="super"
export primer="OPR"
export sample="Zmock"
export minlength="1"
export sample_size="sample_size_all_reads"
export file_id="all_samples_subsampled_all_reads.fastq"

export TRUNC_LEN=0
export THREADS=16


#############################################
# FIXED PROJECT PATHS
#############################################

Base_dir="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon"
Pipeline_dir="${Base_dir}/6_snakemake_usearch_run_all"
Output_root="${Pipeline_dir}/DADA2_ASV_generation"


source ./2_scripts/qiime2_ASVs_generation.sh
