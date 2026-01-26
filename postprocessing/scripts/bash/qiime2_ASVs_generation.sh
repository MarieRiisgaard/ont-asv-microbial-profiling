#!/usr/bin/env bash
set -euo pipefail

##############################################
## USER INPUTS — CHANGE ONLY THESE
##############################################
#
#basecalling="pacbio"              # pacbio | super | hyper
#primer="V4"                       # V4 | V1V3 | V1V8 | OPR
#sample="Zmock"                    # Zmock | ZymoMock | etc.
#minlength="1"                     # e.g. "1"
#sample_size="sample_size_1000"    # sample_size_1000 | sample_size_all_reads
#file_id="all_samples_subsampled_1000.fastq"
#
## Optional trimming
#TRUNC_LEN=0
#THREADS=16
#
##############################################
## FIXED PROJECT PATHS
##############################################
#
#Base_dir="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon"
#Pipeline_dir="${Base_dir}/6_snakemake_usearch_run_all"
#Output_root="${Pipeline_dir}/DADA2_ASV_generation"

#############################################
# BUILD INPUT PATH AUTOMATICALLY
#############################################

# Folder structure example:
# pacbio_all_V4_Zmock_minlength_1/output/03-subsample/sample_size_1000/*.fastq
run_folder="${basecalling}_all_${primer}_${sample}_minlength_${minlength}"

INPUT_FASTQ="${Pipeline_dir}/${run_folder}/output/03-subsample/${sample_size}/${file_id}"

echo "=========================================="
echo " Assembled FASTQ path:"
echo "  $INPUT_FASTQ"
echo "=========================================="

if [[ ! -f "$INPUT_FASTQ" ]]; then
    echo "❌ ERROR: FASTQ file not found!"
    exit 1
fi

#############################################
# BUILD OUTPUT DIRECTORY AUTOMATICALLY
#############################################

OUTDIR="${Output_root}/${basecalling}/${primer}/${sample_size}"
mkdir -p "$OUTDIR"

echo "Output directory:"
echo "  $OUTDIR"
echo "=========================================="

#############################################
# ACTIVATE QIIME2 ENV
#############################################
source "/usr/local/Miniforge3-25.3.1-0-Linux-x86_64/etc/profile.d/conda.sh"
conda activate qiime2-amplicon-2025.10

#############################################
# CREATE MANIFEST
#############################################

MANIFEST="${OUTDIR}/manifest.csv"
echo "sample-id,absolute-filepath,direction" > "$MANIFEST"
echo "sample1,$(readlink -f "$INPUT_FASTQ"),forward" >> "$MANIFEST"

#############################################
# IMPORT INTO QIIME2
#############################################

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path "$MANIFEST" \
  --output-path "${OUTDIR}/demux.qza" \
  --input-format SingleEndFastqManifestPhred33

#############################################
# RUN DADA2
#############################################

qiime dada2 denoise-single \
  --i-demultiplexed-seqs "${OUTDIR}/demux.qza" \
  --p-trim-left 0 \
  --p-trunc-len "$TRUNC_LEN" \
  --p-n-threads "$THREADS" \
  --o-representative-sequences "${OUTDIR}/repseqs.qza" \
  --o-table "${OUTDIR}/asv_table.qza" \
  --o-denoising-stats "${OUTDIR}/stats.qza"\
  --o-base-transition-stats "${OUTDIR}/base_transition_stats.qza"

#############################################
# EXPORT RESULTS
#############################################

# ASV seqs
qiime tools export \
  --input-path "${OUTDIR}/repseqs.qza" \
  --output-path "${OUTDIR}/repseqs"

# ASV table biom
qiime tools export \
  --input-path "${OUTDIR}/asv_table.qza" \
  --output-path "${OUTDIR}/table"

# Convert biom -> TSV
biom convert \
  -i "${OUTDIR}/table/feature-table.biom" \
  -o "${OUTDIR}/table/feature-table.tsv" \
  --to-tsv

#############################################
echo "DONE!"
echo "ASVs saved in:"
echo "  ${OUTDIR}/repseqs/dna-sequences.fasta"
echo "Counts in:"
echo "  ${OUTDIR}/table/feature-table.tsv"
echo "=========================================="



# -------------------------------------------------------------
#  Classify ASVs by mapping to Zymo reference amplicons
# -------------------------------------------------------------

base="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon"
ref="${base}/1_data/zymo_communities/Zymo_Ziels/reference_amplicons/${primer}_bacterial_amplicons.fasta"

asv_fa="${OUTDIR}/repseqs/dna-sequences.fasta"

conda activate snakemake_usearch


# --- Exact matches (100% identity) ---
usearch -search_exact "${asv_fa}" \
    -db "${ref}" \
    -strand both \
    -userout "${OUTDIR}/ASV_exact.tsv" \
    -userfields query+target+id \
    -quiet

# --- Weak / approximate matches (>= 20% identity) ---
usearch -usearch_global "${asv_fa}" \
    -db "${ref}" \
    -id 0.2 \
    -strand both \
    -userout "${OUTDIR}/ASV_weak.tsv" \
    -userfields query+target+id+qcov+tcov+diffs \
    -quiet \
    -maxaccepts 0 -maxrejects 0