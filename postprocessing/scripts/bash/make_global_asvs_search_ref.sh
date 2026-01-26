

#!/usr/bin/env bash
set -euo pipefail


cd /home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/6_snakemake_usearch_run_all
#conda activate snakemake_usearch

### CHANGE MANUALLY ###
primer="V1V8"
include=("pacbio" "super" "hyper")   # Bash array (the folder name include either hyper, pacbio, or super)
########################
outdir="global_ASVs_Ziels/${include[*]}/global_ASVs_${primer}"
#outdir="global_ASVs_updated/${include}/global_ASVs_${primer}"
mkdir -p "$outdir"

echo "▶ Primer = $primer"


# -------------------------------------------------------------
#  Find zOTU files for selected primer and include methods
# -------------------------------------------------------------
zotu_files_tmp="${outdir}/file_list.txt"

# Build list of directories that match both primer and include keyword
zm_dirs=()
for inc in "${include[@]}"; do
    while IFS= read -r dir; do
        zm_dirs+=("$dir")
    done < <(find . -maxdepth 1 -type d -name "*${inc}*${primer}*")
done

# Print what we found
echo "▶ Found directories:"
printf '  %s\n' "${zm_dirs[@]}"

# Collect all zOTU files
> "$zotu_files_tmp"
for d in "${zm_dirs[@]}"; do
    find "$d" -type f -name "zOTUs.fa.filtered" >> "$zotu_files_tmp"
done

echo "▶ zOTU files list saved to: $zotu_files_tmp"

# -------------------------------------------------------------
#  Combine all with annotated headers
# -------------------------------------------------------------
out_annot="${outdir}/all_annotated.fa"
> "$out_annot"

while IFS= read -r f; do
    subset="$(basename "$(dirname "$f")")"  # e.g. sample_size_500
    # Go up 4 levels: file → subset → 04-denoise → output → RUNFOLDER
    runfolder="$(basename "$(dirname "$(dirname "$(dirname "$(dirname "$f")")")")")"

    while IFS= read -r line; do
        if [[ "$line" == ">"* ]]; then
            orig="${line#>}"  # remove '>'
            echo ">${orig}_${runfolder}_${subset}" >> "$out_annot"
        else
            echo "$line" >> "$out_annot"
        fi
    done < "$f"
done < "$zotu_files_tmp"

echo "✅ Combined annotated FASTA written to: $out_annot"



# -------------------------------------------------------------
#  Deduplicate to global ASVs
# -------------------------------------------------------------
usearch -fastx_uniques "${outdir}/all_annotated.fa" \
        -fastaout "${outdir}/global_ASVs.fa" \
        -sizeout -relabel "${primer}_GZOTU"

# -------------------------------------------------------------
#  Map each annotated seq -> global ASV ID
# -------------------------------------------------------------
usearch -search_exact "${outdir}/all_annotated.fa" \
        -db "${outdir}/global_ASVs.fa" \
        -strand plus \
        -userout "${outdir}/map_raw.tsv" \
        -userfields query+target

# Tidy mapping
# query looks like: ORIGID_RUNFOLDER_SAMPLE_SIZE_XXXX
# We want: seq, runfolder, subset, global_asv
awk -F'\t' '
BEGIN{
  OFS="\t";
  print "seq","runfolder","subset","global_asv"
}
{
  q=$1; g=$2;
  # Capture: (anything)_(runfolder)_(sample_size_...)
  # runfolder may contain underscores; subset always starts with sample_size_
  if (match(q, /^(.*)_(.*)_(sample_size_[^[:space:]]+)$/, a)) {
    runfolder=a[2];
    subset=a[3];
  } else {
    runfolder="NA";
    subset="NA";
  }
  print q, runfolder, subset, g
}' "${outdir}/map_raw.tsv" > "${outdir}/map_final.tsv"

echo "✅ Mapping done"
echo "📁 ${outdir}/map_final.tsv"

# -------------------------------------------------------------
#  Taxonomic classification (SINTAX, MiDAS DB)
# -------------------------------------------------------------
usearch -sintax \
    "${outdir}/global_ASVs.fa" \
    -db "/databases/midas/MiDAS5.3_20240320/FLASVs_w_sintax.fa" \
    -tabbedout "${outdir}/global_ASVs_sintax.sintax" \
    -strand both \
    -sintax_cutoff 0.8 \
    -threads 4

# Keep file numerically sorted (stable id order)
sort -V "${outdir}/global_ASVs_sintax.sintax" -o "${outdir}/global_ASVs_sintax.sintax"

echo "✅ SINTAX classification written to:"
echo "   ${outdir}/global_ASVs_sintax.sintax"

# -------------------------------------------------------------
#  Reference matching to Zymo amplicons (exact + weak)
#  (Run on global_ASVs.fa, not per zOTU file)
# -------------------------------------------------------------
base="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon"
ref="${base}/1_data/zymo_communities/Zymo_Ziels/reference_amplicons/${primer}_bacterial_amplicons.fasta"

# --- Exact matches (100% identity) ---
usearch -search_exact "${outdir}/global_ASVs.fa" \
    -db "${ref}" \
    -strand both \
    -userout "${outdir}/global_ASVs_exact.tsv" \
    -userfields query+target+id \
    -quiet

# --- Weak / approximate matches (>= 0.2 identity) ---
usearch -usearch_global "${outdir}/global_ASVs.fa" \
    -db "${ref}" \
    -id 0.2 \
    -strand both \
    -userout "${outdir}/global_ASVs_weak.tsv" \
    -userfields query+target+id+qcov+tcov+diffs \
    -quiet \
    -maxaccepts 0 -maxrejects 0

echo "✅ Reference matching done"
echo "📁 Exact: ${outdir}/global_ASVs_exact.tsv"
echo "📁 Weak : ${outdir}/global_ASVs_weak.tsv"




