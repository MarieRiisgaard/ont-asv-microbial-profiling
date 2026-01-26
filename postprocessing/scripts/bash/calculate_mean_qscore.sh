#!/usr/bin/env bash
conda activate evaluate_ASVs

# ---------------- USER SETTINGS ----------------
base="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/6_snakemake_usearch_run_all"
outdir="${base}/read_quality"
mkdir -p "$outdir"

methods=("pacbio" "super")
primers=("V4" "V1V3" "V1V8" "OPR")
samples=("Zmock" "Zfecal" "WWTP" "AD" "Soil")

# Read-depth rules
reads_V4=10000
reads_long=100000
# ------------------------------------------------

for method in "${methods[@]}"; do
  for primer in "${primers[@]}"; do
    for sample in "${samples[@]}"; do

      # Select read depth based on primer
      if [[ "$primer" == "V4" ]]; then
        reads=$reads_V4
      else
        reads=$reads_long
      fi

      # Construct expected FASTQ path
      fastq="${base}/${method}_all_${primer}_${sample}_minlength_1/output/03-subsample/sample_size_${reads}/all_samples_subsampled_${reads}.fastq"

      # Check if file exists
      if [[ ! -f "$fastq" ]]; then
        echo "[SKIP] Missing FASTQ: $fastq"
        continue
      fi

      # Output filename
      outfile="${outdir}/${method}_${primer}_${sample}_${reads}_Qbins.tsv"

      echo "[RUN] Processing $method $primer $sample ($reads reads)"

      # Extract Q-score bins
      seqkit fx2tab -n -q -l "$fastq" \
        | awk '
            ($4 ~ /^[0-9.]+$/){
                q = int($4)
                counts[q]++
            }
            END {
                print "Q_bin\tcount"
                for (q in counts) print q "\t" counts[q]
            }' \
        | sort -k1,1n \
        > "$outfile"

      echo "[DONE] → $outfile"
      echo ""

    done
  done
done







# ------------------------------------------------

