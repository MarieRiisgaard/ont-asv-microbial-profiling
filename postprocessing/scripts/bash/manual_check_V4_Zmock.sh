
# Only the entries including samonella
V4_Zmock_ref=/home/bio.aau.dk/hd95lp/Projects/RETHINK/ont-asv-microbial-profiling/postprocessing/data/Zmock_in_silico_amplicons/V4_bacterial_amplicons.fasta   # output: V4_Zmock_ref_samonella
#the format: 
>Salmonella_5;size=5
TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTGTCAAGTCGGATGTGAAATCCC
CGGGCTCAACCTGGGAACTGCATTCGAAACTGGCAGGCTTGAGTCTTGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTG
AAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACAAAGACTGACGCTCAGGTGCGAAAGCGTG
GGGAGCAAACAGG
>Enterococcus_4;size=4
TACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTTCTTAAGTCTGATGTGAAAGCCC
CCGGCTCAACCGGGGAGGGTCATTGGAAACTGGGAGACTTGAGTGCAGAAGAGGAGAGTGGAATTCCATGTGTAGCGGTG
AAATGCGTAGATATATGGAGGAACACCAGTGGCGAAGGCGGCTCTCTGGTCTGTAACTGACGCTGAGGCTCGAAAGCGTG
GGGAGCAAACAGG

All output in /home/bio.aau.dk/hd95lp/Projects/RETHINK/ont-asv-microbial-profiling/postprocessing/results/manual_check_V4_zmock

V4_Zmock_ref_samonella=subset to only asamonella (there are two seqs)



subset to 100000 = V4_primer_trimmed_reads_PB subset (for testing)


usearch -search_exact "$V4_primer_trimmed_reads_PB" \
        -db "$V4_primer_trimmed_reads_PB" \
        -strand plus \
        -userout "${outdir}/map_raw.tsv" \
        -userfields query+target

# I want to know 
how many reads maps to each of the reference seq
(the V4_primer_trimmed_reads_PB is merged for four barcode) --> tell me the result for both the merged and each barcode (I want to test this using the subset for computational speed)
I need to run it seperately for the PB and ONT 


conda activate snakemake_usearch
outdir=/home/bio.aau.dk/hd95lp/Projects/RETHINK/ont-asv-microbial-profiling/postprocessing/results/manual_check_V4_zmock
mkdir -p "$outdir"

# Extract the two Salmonella references
V4_Zmock_ref=/home/bio.aau.dk/hd95lp/Projects/RETHINK/ont-asv-microbial-profiling/postprocessing/data/Zmock_in_silico_amplicons/V4_bacterial_amplicons.fasta

awk '/^>/ {keep=(tolower($0) ~ /salmonella/)} keep' \
    "$V4_Zmock_ref" > "$outdir/V4_Zmock_ref_Salmonella.fasta"

V4_primer_trimmed_reads_ONT=/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/6_snakemake_usearch_unoise_minsize_test/super_all_V4_Zmock_unoise_minsize_8/output/03-subsample/sample_size_all_reads/all_samples_subsampled_all_reads.fastq

V4_primer_trimmed_reads_PB=/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/6_snakemake_usearch_unoise_minsize_test/pacbio_all_V4_Zmock_unoise_minsize_8/output/03-subsample/sample_size_all_reads/all_samples_subsampled_all_reads.fastq


# Choose PB or ONT
platform=PB
reads="$V4_primer_trimmed_reads_PB"

# For ONT, use these instead:
# platform=ONT
# reads="$V4_primer_trimmed_reads_ONT"

# Find full-length perfect matches
usearch -search_exact "$V4_primer_trimmed_reads_PB" \
    -db "$outdir/V4_Zmock_ref_Salmonella.fasta" \
    -strand plus \
    -userout "$outdir/${platform}_map_raw.tsv" \
    -userfields query+target

platform=ONT
usearch -search_exact "$V4_primer_trimmed_reads_ONT" \
    -db "$outdir/V4_Zmock_ref_Salmonella.fasta" \
    -strand plus \
    -userout "$outdir/${platform}_map_raw.tsv" \
    -userfields query+target

# Count matches per reference across all barcodes
cut -f2 "$outdir/${platform}_map_raw.tsv" | sort | uniq -c