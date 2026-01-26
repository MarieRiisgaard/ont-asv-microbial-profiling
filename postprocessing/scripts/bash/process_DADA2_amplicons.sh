#################################
BASE="/home/bio.aau.dk/hd95lp/Projects/RETHINK/ASV_resolved_community_with_ONT/postprocessing/results/DADA2_ASV_generation"
OUTDIR="${BASE}/global_ASVs_DADA2_"
OUTBASE="${BASE}/global_ASVs_DADA2_"

mkdir -p "$OUTDIR"

PRIMERS=("V4" "V1V3" "V1V8" "OPR")

for PRIMER in "${PRIMERS[@]}"; do
  OUT="${OUTDIR}/merged_${PRIMER}_dna-sequences.fasta"
  : > "$OUT"

  find "$BASE" -type f -path "*/${PRIMER}/*/repseqs/dna-sequences.fasta" -print0 \
  | sort -z \
  | while IFS= read -r -d '' f; do

      # Make path relative to BASE, then split into components
      rel="${f#${BASE}/}"   # e.g. pacbio/OPR/sample_size_100/repseqs/dna-sequences.fasta
      IFS='/' read -r platform primerset sample_size repdir _ <<< "$rel"

      prefix="${platform}_${primerset}_${sample_size}"

      awk -v pfx="$prefix" '
        BEGIN { OFS=""; c=0; }
        /^>/ { c++; print ">seq", c, "_", pfx; next }
        { gsub(/[[:space:]]+/, "", $0); if (length($0)>0) print $0 }
      ' "$f" >> "$OUT"
    done

  echo "Wrote: $OUT"
done



#######################

# Your merged files from earlier step:
# ${OUTBASE}/merged_V4_dna-sequences.fasta
# ${OUTBASE}/merged_V1V3_dna-sequences.fasta
# ${OUTBASE}/merged_V1V8_dna-sequences.fasta
# ${OUTBASE}/merged_OPR_dna-sequences.fasta

PRIMERS=("V4" "V1V3" "V1V8" "OPR")

# Reference directory for in silico amplicons
REFDIR="/home/bio.aau.dk/hd95lp/Projects/RETHINK/ASV_resolved_community_with_ONT/postprocessing/data/Zmock_in_silico_amplicons"

# Chimera-check output base
UCHIME_BASE="${BASE}/uchime_ref"
mkdir -p "$UCHIME_BASE"

# If you want to run inside conda:
# conda activate snakemake_usearch

for primer in "${PRIMERS[@]}"; do
  echo "======================================"
  echo "Primer: ${primer}"
  echo "======================================"

  merged="${OUTBASE}/merged_${primer}_dna-sequences.fasta"
  outdir="${OUTBASE}/global_ASVs_${primer}"
  mkdir -p "$outdir"

  if [[ ! -s "$merged" ]]; then
    echo "WARNING: Missing or empty merged file: $merged"
    echo "Skipping ${primer}."
    continue
  fi

  # -------------------------------------------------------------
  # Deduplicate to global ASVs
  # -------------------------------------------------------------
  global_asvs="${outdir}/global_ASVs_${primer}.fa"

  usearch -fastx_uniques "$merged" \
          -fastaout "$global_asvs" \
          -sizeout \
          -relabel "${primer}_GZOTU"

  # -------------------------------------------------------------
  # Map each annotated seq -> global ASV ID
  # -------------------------------------------------------------
  if [[ -s "$merged" ]]; then
    usearch -search_exact "$merged" \
            -db "$global_asvs" \
            -strand plus \
            -userout "${outdir}/map_raw.tsv" \
            -userfields query+target
  else
    echo "NOTE: all_annotated.fa not found at: $merged"
    echo "Skipping mapping step for ${primer}."
  fi

  # -------------------------------------------------------------
  # Exact matches vs in silico reference for this primer
  # -------------------------------------------------------------
  ref="${REFDIR}/${primer}_bacterial_amplicons.fasta"

  if [[ -s "$ref" ]]; then
    usearch -search_exact "$global_asvs" \
            -db "$ref" \
            -strand both \
            -userout "${outdir}/global_ASVs_exact.tsv" \
            -userfields query+target+id \
            -quiet
  else
    echo "NOTE: Reference not found for ${primer}: $ref"
    echo "Skipping exact-match step for ${primer}."
  fi

  # -------------------------------------------------------------
  # Make CSV: global OTU ID, Primer, sequence length
  # (parses FASTA, removes any whitespace, counts sequence chars)
  # -------------------------------------------------------------
  csv="${outdir}/global_ASVs_${primer}_lengths.csv"
  awk -v primer="$primer" '
    BEGIN { OFS=","; print "global_otu_id,primer,sequence_length"; }
    /^>/ {
      if (id != "") { print id, primer, seqlen; }
      id = substr($0,2);
      seqlen = 0;
      next
    }
    {
      gsub(/[[:space:]]+/, "", $0);
      seqlen += length($0);
    }
    END {
      if (id != "") { print id, primer, seqlen; }
    }
  ' "$global_asvs" > "$csv"

  # -------------------------------------------------------------
  # Chimera check (uchime2_ref)
  # -------------------------------------------------------------
  # Output paths
  chimdir="${UCHIME_BASE}/${primer}"
  mkdir -p "$chimdir"

  if [[ -s "$ref" ]]; then
    usearch -uchime2_ref "$global_asvs" \
      -db "$ref" \
      -strand plus \
      -mode balanced \
      -chimeras "${chimdir}/asvs_chimeras_ref.fa" \
      -uchimeout "${chimdir}/asvs_uchime_ref.txt"
  else
    echo "NOTE: Reference not found for chimera check for ${primer}: $ref"
    echo "Skipping chimera check for ${primer}."
  fi

  echo "Done: ${primer}"
  echo "Outputs in: ${outdir}"
done



############
############