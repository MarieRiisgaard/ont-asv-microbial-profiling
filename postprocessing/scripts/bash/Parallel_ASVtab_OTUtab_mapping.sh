#!/usr/bin/bash -l
#SBATCH --job-name=MFD_fSSU_vs_MFD_ssu_database_v1.0_ASVtab_mapping
#SBATCH --output=6_snakemake_usearch_run_all/OTU_tab_global_ASVs/slurm/job_%j_%x.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=TIME_LIMIT_80
#SBATCH --ntasks-per-node=1
#SBATCH --partition=general
#SBATCH --cpus-per-task=40
#SBATCH --mem=100G
#SBATCH --time=1-22:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mari@bio.aau.dk

#----- THIS script is copied from MD
# /home/bio.aau.dk/md/Common_scripts/Amplicons/Parallel_ASVtab_exact_mapping.sh


# --------- I have added ---------------


# Exit on first error and if any variables are unset
set -eu
max_threads="$(nproc)"

otutab_parallel() {
  set -eu

  local OPTIND input="" zotus="" output="" threads=8 splits=10
  while getopts ":i:z:o:t:s:" opt; do
    case $opt in
      i ) input=$OPTARG ;;
      z ) zotus=$OPTARG ;;
      o ) output=$OPTARG ;;
      t ) threads=$OPTARG ;;
      s ) splits=$OPTARG ;;
    esac
  done

  if [[ -z "$input" || -z "$zotus" || -z "$output" ]]; then
      echo "Usage: otutab_parallel -i reads.fa -z zotus.fa -o otutab.txt -t threads -s splits"
      exit 1
  fi

  tmpdir="${output}_tmp"
  rm -rf "$tmpdir"
  mkdir -p "$tmpdir"

  echo ">>> Splitting input into $splits chunks..."
  usearch -fastx_split "$input" \
      -splits "$splits" \
      -outname "$tmpdir/chunk_@.fa" \
      -quiet

  job_threads=$(( threads / splits ))
  (( job_threads < 1 )) && job_threads=1

  echo ">>> Running $splits parallel OTU-tab jobs..."
  find "$tmpdir" -name "chunk_*.fa" | \
      parallel --halt soon,fail=1 \
          usearch -otutab {} \
          -zotus "$zotus" \
          -otutabout {.}.otu.txt \
          -threads "$job_threads" \
          -quiet

  echo ">>> Merging OTU tables using Python (safe method)..."

python3 <<EOF
import pandas as pd
import glob

tables = glob.glob("${tmpdir}/*.otu.txt")

dfs = []
for tbl in tables:
    # Important: KEEP the '#OTU ID' header line
    df = pd.read_csv(tbl, sep="\t", header=0, comment=None)

    # Rename '#OTU ID' → 'OTU'
    first_col = df.columns[0]
    df = df.rename(columns={first_col: "OTU"})

    # Set OTU as index
    df = df.set_index("OTU")

    dfs.append(df)

# Sum all tables element-wise
merged = dfs[0].copy()
for df in dfs[1:]:
    merged = merged.add(df, fill_value=0)

# Sort for consistent output
merged = merged.sort_index(axis=0).sort_index(axis=1)

# Save final table
merged.to_csv("${output}", sep="\t")
EOF


  echo ">>> Merged OTU table written to: ${output}"
  # rm -rf "$tmpdir"
}


###########################################

source "/usr/local/Miniforge3-25.3.1-0-Linux-x86_64/etc/profile.d/conda.sh"

########## USER VARIABLES ##########
method="super"        # super or pacbio
primer="OPR"           # V4, V1V3, V1V8, OPR
sample="Zfecal"        # dataset
reads="all_reads"     # can be e.g. 50000, all_reads etc.
####################################

base="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/6_snakemake_usearch_run_all"
input_fastq="${base}/${method}_all_${primer}_${sample}_minlength_1/output/03-subsample/sample_size_${reads}/all_samples_subsampled_${reads}.fastq"
outdir="${base}/OTU_tab_global_ASVs/${reads}/${sample}/${primer}"
mkdir -p "$outdir"

output_fasta="${outdir}/${method}_${sample}_${primer}_all_samples_subsampled_reads.fa"
tmp_fasta="${output_fasta%.fa}_raw.fa"


ASVs_fasta="${base}/global_ASVs_Ziels/pacbio_super/global_ASVs_${primer}/global_ASVs.fa"
csv_ids="${base}/global_ASVs_Ziels/pacbio_super/${primer}_${sample}_all_reads.csv"
#output_fa="${base}/global_ASVs_Ziels/pacbio_super/global_ASVs_${primer}_${sample}/selected_ASVs_all_reads.fa"
outdir_asv="${base}/global_ASVs_Ziels/pacbio_super/global_ASVs_${primer}"
mkdir -p "$outdir_asv"

output_fa="${base}/global_ASVs_Ziels/pacbio_super/global_ASVs_${primer}/selected_ASVs_${sample}_all_reads.fa"


echo "Extracting sequences from:"
echo "  FASTA: $ASVs_fasta"
echo "  CSV:   $csv_ids"
echo "  OUT:   $output_fa"

# Add semicolon to each ID to mimic FASTA header prefix
sed 's/$/;/' "$csv_ids" > tmp_headers.txt

# Extract sequences whose header (before first ;) matches any ID
awk -v hdrs="tmp_headers.txt" '
  BEGIN {
    # Load IDs into lookup table
    while (getline < hdrs) {
      wanted[$0] = 1
    }
    RS=">"; ORS=""
  }

  NR == 1 { next }  # skip initial empty split

  {
    # header is everything before the first semicolon
    split($1, a, ";")
    header = a[1] ";"  # reconstruct the part we want to match

    if (header in wanted)
      print ">" $0
  }
' "$ASVs_fasta" > "$output_fa"

echo "Extracted $(grep -c "^>" "$output_fa") sequences into:"
echo "  $output_fa"

# Build exact-match regex: V4_GZOTU3; | V4_GZOTU20; | ...
#pattern=$(sed 's/$/;/' "$csv_ids" | paste -sd '|' -)
#
## Extract exact matches into output_fa
#awk -v pat="$pattern" '
#  BEGIN {RS=">"; ORS=""}
#  NR == 1 {next}
#
#  {
#    # header up to first semicolon
#    split($1, h, ";")
#    header = h[1] ";"
#
#    if (header ~ ("^(" pat ")$"))
#      print ">" $0
#  }
#' "$ASVs_fasta" > "$output_fa"
#
#echo "Done. Extracted $(grep -c '^>' "$output_fa") sequences into:"
##echo "  $output_fa"






# Threads
max_threads=$(nproc --ignore=1)

echo ">>> Running FASTQ→FASTA for $method - $primer - $sample"

## Convert FASTQ to FASTA
conda activate porechop_env

# 1) Convert FASTQ → FASTA
seqtk seq -a "$input_fastq" > "$tmp_fasta"
# 2) Keep only >barcodeXX (remove _something and trailing text)
awk '
/^>/ {
    header=$1
    sub(/_.*/, "", header)   # remove underscore and after
    sub(/ .*/, "", header)   # remove space and after
    print header
    next
}
{ print }
' "$tmp_fasta" > "$output_fasta"
# 3) Remove temp file
rm "$tmp_fasta"

conda activate snakemake_usearch

# ---- RUN PARALLEL OTU-TAB CREATION ---- #

echo ">>> Running parallel OTU-tab: $method - $primer - $sample"

otutab_parallel \
    -i "$output_fasta" \
    -z "$output_fa" \
    -o "${outdir}/${method}_OTUtab.txt" \
    -t "$max_threads" \
    -s 10

echo "✅ Finished OTU-tab: ${method} - ${primer} - ${sample}"
echo "Output table: ${outdir}/${method}_OTUtab.txt"


# sbatch 2_scripts/Parallel_ASVtab_OTUtab_mapping.sh

#usearch -tabbed2otutab temp.m8 -output /home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/6_snakemake_usearch_run_all/ASV_tab_exact_matches/Zmock/V4/super_ASVtab.txt