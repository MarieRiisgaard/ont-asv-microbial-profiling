#!/bin/bash
set -euo pipefail

#conda activate evaluate_ASVs

# -----------------------------
# User-defined paths
# -----------------------------
#full_16S_dir="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/1_data/zymo_communities/Zfecal/RNA"
full_16S_dir="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/1_data/zymo_communities/Zymo_Ziels"
#full_16S_dir="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/1_data/zymo_communities/Zymo_Ziels/zymo-ref-uniq_2019-10-28.fa"

#output_dir="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/1_data/zymo_communities/Zfecal/amplicons"
output_dir="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/1_data/zymo_communities/Zymo_Ziels/amplicons"


mkdir -p "$output_dir"

# -----------------------------
# Primer sets to process
# -----------------------------
for primers in V1V8 V4 V1V3; do
    echo "🔹 Processing primer set: $primers"

    case $primers in
        V1V8)
            FWD="AGRGTTYGATYMTGGCTCAG"
            REV="TGYACWCACCGCCCGTC"
            ;;
        V4)
            FWD="GTGYCAGCMGCCGCGGTAA"
            REV="ATTAGAWACCCBDGNTAGTCC"
            ;;
        V1V3)
            FWD="AGRGTTYGATYMTGGCTCAG"
            REV="CCAGCAGCCGCGGTAAT"
            ;;
        *)
            echo "❌ Unknown primer set: $primers"
            exit 1
            ;;
    esac

    # Output paths
    trimmed_out="${output_dir}/primer_${primers}_bacterial_amplicons.fasta"
    untrimmed_out="${output_dir}/primer_${primers}_untrimmed.fasta"
    log_file="${output_dir}/primer_${primers}_cutadapt_log.txt"

    # Combine all 16S sequences into one temporary file
    tmp_input="${output_dir}/tmp_all_16S_${primers}.fasta"
    #cat "$full_16S_dir"/*_dereplicated_renamed.fasta > "$tmp_input"
    cat "$full_16S_dir"/*.fa > "$tmp_input"


    echo "Running cutadapt for $primers ..."
    cutadapt \
        -g "${FWD}...${REV}" \
        -o "$trimmed_out" \
        --revcomp \
        --untrimmed-output "$untrimmed_out" \
        "$tmp_input" \
        -j 4 > "$log_file" 2>&1

    echo "✅ Primer trimming finished for $primers. Log: $log_file"

    # Dereplicate
    output_derep="${output_dir}/primer_${primers}_bacterial_amplicons_dereplicated.fasta"
    vsearch --derep_fulllength "$trimmed_out" \
            --output "$output_derep" \
            --sizeout \
            --threads 4

    echo "✅ Dereplication finished for $primers -> $output_derep"
    rm "$tmp_input"

    echo ""
done

echo "🎉 All primer sets processed successfully!"


# -------------------------------------------------------------
#  Manually make the OPR reference amplicons (Zymo) 
# -------------------------------------------------------------

cd "/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/1_data/zymo_communities/Zymo_Ziels/make_OPR_cuts" \
  || { echo "❌ Folder not found"; exit 1; }

# Input FASTA of full-length Zymo genomes' 16S sequences
OPR_amplicons_zymo="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/1_data/zymo_communities/Zymo_Ziels/zymo-ref-uniq_2019-10-28.fa"

# -------------------------------------------------------------
#  Step 1 — Rename sequences with species names
# -------------------------------------------------------------
python3 <<EOF
#!/usr/bin/env python3
import re
from collections import defaultdict

input_fasta = "$OPR_amplicons_zymo"
output_fasta = input_fasta.replace(".fasta", "_renamed.fasta")

mapping = {
    "CP039755.1": "Bacillus_subtilis",
    "CP039752.1": "Enterococcus_faecalis",
    "CP039753.1": "Escherichia_coli",
    "CP039750.1": "Limosilactobacillus_fermentum",
    "CP039751.1": "Listeria_monocytogenes",
    "CP039749.1": "Pseudomonas_aeruginosa",
    "CP039757.1": "Salmonella_enterica",
    "CP039759.1": "Staphylococcus_aureus",
}

counters = defaultdict(int)

def clean_header(header):
    acc_match = re.match(r">(CP\d+\.\d+)", header)
    if not acc_match:
        return header
    acc = acc_match.group(1)
    species = mapping.get(acc, acc)
    counters[species] += 1
    idx = counters[species]
    return f">{species}_{idx}"

with open(input_fasta) as fin, open(output_fasta, "w") as fout:
    header = None
    seq_lines = []
    for line in fin:
        line = line.strip()
        if line.startswith(">"):
            if header:
                fout.write(header + "\n" + "".join(seq_lines) + "\n")
            header = clean_header(line)
            seq_lines = []
        else:
            seq_lines.append(line)
    if header:
        fout.write(header + "\n" + "".join(seq_lines) + "\n")

print(f"✅ Renamed FASTA written to: {output_fasta}")
EOF

# Update variable to point to renamed file
#OPR_amplicons_renamed="${OPR_amplicons_zymo/.fasta/_renamed.fasta}"
OPR_amplicons_renamed="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/1_data/zymo_communities/Zymo_Ziels/zymo-ref-uniq_2019-10-28.fa"
# -------------------------------------------------------------
#  Step 2 — Extract OPR region with cutadapt
# -------------------------------------------------------------
FWD="AGRGTTYGATYMTGGCTCAG"
REV="GTTTGGCACCTCGATGTCG"

cutadapt \
  -g "${FWD}...${REV}" \
  -o "trimmed_OPR_Zmock.fasta" \
  --revcomp \
  --untrimmed-output "untrimmed_OPR_Zmock.fasta" \
  "$OPR_amplicons_renamed" \
  -j 4 > "trim_log.txt" 2>&1

# -------------------------------------------------------------
#  Step 3 — Dereplicate trimmed amplicons
# -------------------------------------------------------------
output_derep="OPR_bacterial_amplicons_dereplicated.fasta"
vsearch --derep_fulllength "trimmed_OPR_Zmock.fasta" \
        --output "$output_derep" \
        --sizeout \
        --fasta_width 0 \
        --threads 4

echo "✅ Done — dereplicated output: $output_derep"
