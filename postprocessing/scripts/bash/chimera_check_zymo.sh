

conda activate snakemake_usearch

zymo_ASVs="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/6_snakemake_usearch_run_all/global_ASVs_Ziels/pacbio_super/global_ASVs_OPR/selected_ASVs_all_reads.fa"
ref_seq="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/1_data/zymo_communities/Zymo_Ziels/reference_amplicons/OPR_bacterial_amplicons.fasta"
outdir="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/6_snakemake_usearch_run_all/uchime_ref/"

# Create output directory
mkdir -p "$outdir/OPR"

usearch -uchime2_ref "$zymo_ASVs" \
  -db "$ref_seq" \
  -strand plus \
  --mode balanced \
  -chimeras "$outdir/OPR/asvs_chimeras_ref.fa" \
  -uchimeout "$outdir/OPR/asvs_uchime_ref.txt"

zymo_ASVs="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/6_snakemake_usearch_run_all/global_ASVs_Ziels/pacbio_super/global_ASVs_V4/selected_ASVs_all_reads.fa"
ref_seq="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/1_data/zymo_communities/Zymo_Ziels/reference_amplicons/V4_bacterial_amplicons.fasta"
# Create output directory
mkdir -p "$outdir/V4"


usearch -uchime2_ref "$zymo_ASVs" \
  -db "$ref_seq" \
  -strand plus \
  -mode balanced \
  -chimeras "$outdir/V4/asvs_chimeras_ref.fa" \
  -uchimeout "$outdir/V4/asvs_uchime_ref.txt"


zymo_ASVs="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/6_snakemake_usearch_run_all/global_ASVs_Ziels/pacbio_super/global_ASVs_V1V3/selected_ASVs_all_reads.fa"
ref_seq="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/1_data/zymo_communities/Zymo_Ziels/reference_amplicons/V1V3_bacterial_amplicons.fasta"
# Create output directory
mkdir -p "$outdir/V1V3"


usearch -uchime2_ref "$zymo_ASVs" \
  -db "$ref_seq" \
  -strand plus \
  -mode balanced \
  -chimeras "$outdir/V1V3/asvs_chimeras_ref.fa" \
  -uchimeout "$outdir/V1V3/asvs_uchime_ref.txt"



zymo_ASVs="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/6_snakemake_usearch_run_all/global_ASVs_Ziels/pacbio_super/global_ASVs_V1V8/selected_ASVs_all_reads.fa"
ref_seq="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/1_data/zymo_communities/Zymo_Ziels/reference_amplicons/V1V8_bacterial_amplicons.fasta"
# Create output directory
mkdir -p "$outdir/V1V8"


usearch -uchime2_ref "$zymo_ASVs" \
  -db "$ref_seq" \
  -strand plus \
  -mode balanced \
  -chimeras "$outdir/V1V8/asvs_chimeras_ref.fa" \
  -uchimeout "$outdir/V1V8/asvs_uchime_ref.txt"
