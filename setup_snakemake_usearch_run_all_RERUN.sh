#!/bin/bash
set -euo pipefail


# -------------------------------------------------------------
# Global settings
# -------------------------------------------------------------
base_dir="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/6_snakemake_usearch_unoise_minsize_test"
#subsample_sizes="[500, 1000, 5000, 10000, 50000, 100000, 200000, 400000, 600000, 800000, 1000000]"
#subsample_sizes="[10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000]" # Mortens suggestion
subsample_sizes="[100, 1000, 10000, 100000, 1000000]" # Mortens suggestion
setting="unoise_minsize_4"

CONDA_ACTIVATE="/usr/local/Miniforge3-25.3.1-0-Linux-x86_64/etc/profile.d/conda.sh"
ENV_NAME="snakemake_usearch"

declare -A datasets=(
    ["super_all"]="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/1_data/primer_test/super_basecalling/all/fastq_barcode_split"
    ["pacbio_all"]="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/1_data/primer_test/pacbio/fastq_barcode_split"
    #["hyper_all"]="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/1_data/primer_test/hyper_basecalling/20251021_hyp_full_demuxed/fastq_barcode_split"
)

# -------------------------------------------------------------
# Helper function to update config keys
# -------------------------------------------------------------
update_or_add () {
    local key=$1
    local value=$2
    local file=config/config.yaml
    if grep -q "^$key:" "$file"; then
        sed -i "s|^$key:.*|$key: $value|" "$file"
    else
        echo "$key: $value" >> "$file"
    fi
}

# -------------------------------------------------------------
# Function to deploy + submit one combination
# -------------------------------------------------------------
run_one () {
    local dataset=$1 primer=$2 sample=$3
    local dataset_root="${datasets[$dataset]}"
    local data_path="$dataset_root/$primer/$sample"
    [ -d "$data_path" ] || { echo "⚠️ Skipping missing: $data_path"; return; }

    local run_id="${dataset}_${primer}_${sample}_${setting}"
    local run_folder="$base_dir/$run_id"
    echo "🚀 Starting $run_id"

    # --- Deploy workflow if needed ---
    if [ ! -f "$run_folder/workflow/Snakefile" ]; then
        snakedeploy deploy-workflow \
            https://github.com/MarieRiisgaard/snakemake_usearch_primertrim_and_subset \
            --target-dir "$run_folder" \
            --branch main \
            --force
    fi

    cd "$run_folder"

    # --- Primer setup ---
    case $primer in
        V1V8) FWD="AGRGTTYGATYMTGGCTCAG"; REV="TGYACWCACCGCCCGTC";;
        V4)   FWD="GTGYCAGCMGCCGCGGTAA"; REV="ATTAGAWACCCBDGNTAGTCC";;
        OPR)  FWD="AGRGTTYGATYMTGGCTCAG"; REV="GTTTGGCACCTCGATGTCG";;
        V1V3) FWD="AGRGTTYGATYMTGGCTCAG"; REV="CCAGCAGCCGCGGTAAT";;
        *) echo "❌ Unknown primer $primer"; return;;
    esac

    # --- Reference file selection ---
    REF=""
    if [ "$sample" = "Zfecal" ]; then
        case $primer in
            V1V8) REF="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/1_data/zymo_communities/Zfecal/amplicons/primer_V1V8_zfecal_bacterial_amplicons_dereplicated.fasta";;
            V4)   REF="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/1_data/zymo_communities/Zfecal/amplicons/primer_V4_zfecal_bacterial_amplicons_dereplicated.fasta";;
            V1V3) REF="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/1_data/zymo_communities/Zfecal/amplicons/primer_V1V3_zfecal_bacterial_amplicons_dereplicated.fasta";;
            OPR)  REF="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/1_data/zymo_communities/Zfecal/genomes/OPR_amplicons_all_dereplicated.fasta";;
        esac
    elif [ "$sample" = "Zmock" ]; then
        REF="/home/bio.aau.dk/hd95lp/Projects/RETHINK/p_nanopore_amplicon/1_data/zymo_communities/Zmock/reference_amplicons/${primer}_bacterial_amplicons.fasta"
    fi

    echo "✅ Reference: ${REF:-none}"

    # --- Update config.yaml ---
    update_or_add "input_dir" "\"$data_path\""
    update_or_add "primers" "\"$FWD...$REV\""
    update_or_add "output_dir" "\"output\""
    update_or_add "filtlong_args" "\"--min_length 1\""
    update_or_add "subsample_sizes" "$subsample_sizes"
    update_or_add "reference_amplicons" "\"$REF\""
    update_or_add "unoise_minsize" "\4\"


    # --- Write SLURM script ---
    mkdir -p logs
    slurm_script="run_snakemake.slurm"
    cat > "$slurm_script" <<EOL
#!/bin/bash
#SBATCH --job-name=${run_id}
#SBATCH --output=${run_folder}/job_%j_%x.out
#SBATCH --cpus-per-task=2
#SBATCH --time=0-12:00:00
#SBATCH --mem=2G
#SBATCH --account=phn

source ${CONDA_ACTIVATE}
conda activate ${ENV_NAME}

snakemake --unlock
snakemake --executor slurm --jobs 96 --keep-going --rerun-incomplete --printshellcmds --set-resources abund_table:runtime=7200

EOL
    chmod +x "$slurm_script"

    echo "🎯 Submitting $run_id"
    sbatch "$slurm_script"

    cd - >/dev/null
}

export -f run_one update_or_add
export base_dir setting subsample_sizes CONDA_ACTIVATE ENV_NAME datasets

# -------------------------------------------------------------
# Parallel launch
# -------------------------------------------------------------
max_parallel=8  # adjust depending on cluster quota

for dataset in "${!datasets[@]}"; do
    for primer in V1V8 V1V3 V4 OPR; do
        for sample in Zmock Zfecal Soil AD WWTP; do
            run_one "$dataset" "$primer" "$sample" &
            # Limit to $max_parallel concurrent background jobs
            while (( $(jobs -r | wc -l) >= max_parallel )); do sleep 10; done
        done
    done
done

wait
echo "✅ All runs submitted successfully!"
