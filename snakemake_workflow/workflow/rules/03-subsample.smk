###############################################################
# 03-subsample.smk — dynamic subsampling after merging all reads
###############################################################

import os, subprocess

# -------------------------------------------------------------
# 🧠 Helper function: dynamically determine valid subsample sizes
# -------------------------------------------------------------
def get_dynamic_subsample_sizes(wildcards):
    """
    Reads the number of reads in the merged FASTQ and returns
    a list of subsample sizes that are smaller than that count,
    plus 'all_reads' for the full dataset.
    """
    merged_fastq = os.path.join(config["tmp_dir"], "merged", "all_samples_trimmed.fastq")

    # 👇 Key line: skip grep if file doesn’t exist yet
    if not os.path.exists(merged_fastq):
        print(f"[INFO] {merged_fastq} not found yet (DAG build). Returning default sizes.")
        return config.get("subsample_sizes", [10000, 20000, 50000, 100000, 200000])

    try:
        nreads = int(subprocess.check_output(["grep", "-c", "^+$", merged_fastq]).strip())
    except Exception:
        nreads = 0

    default_sizes = config.get("subsample_sizes", [10000, 20000, 50000, 100000, 200000])
    valid_sizes = [s for s in default_sizes if s < nreads]

    if nreads > 0:
        valid_sizes.append("all_reads")

    print(f"[INFO] Total reads: {nreads}. Valid subsample sizes: {valid_sizes}")
    return valid_sizes


# -------------------------------------------------------------
# 1️⃣ Concatenate all trimmed reads (across all barcodes/samples)
# -------------------------------------------------------------
rule concat_all_trimmed:
    input:
        expand(
            os.path.join(
                config["tmp_dir"],
                "02-trim_primers",
                "{sample}",
                "{sample}_trimmed.fastq"
            ),
            sample=sample_dirs,
        )
    output:
        os.path.join(config["tmp_dir"], "merged", "all_samples_trimmed.fastq")
    log:
        os.path.join(config["log_dir"], "03-subsample", "concat_all_trimmed.log")
    message:
        "Concatenating all trimmed FASTQs into a single merged file"
    container:
        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
    conda:
        "../envs/snakemake_usearch.yml"
    resources:
        mem_mb=32000,     # 👈 give it 32 GB
        runtime=120,      # minutes or use SLURM default if not handled
    threads: 1
    shell:
        r"""
        exec &> "{log}"
        set -euxo pipefail

        mkdir -p "$(dirname {output})"
        cat {input} > "{output}"

        # Basic QC
        total_reads=$(grep -c '^+$' "{output}" || true)
        echo "Total merged reads: $total_reads"
        if [ "$total_reads" -eq 0 ]; then
            echo "❌ No reads found in merged file!"
            exit 1
        fi
        """


# -------------------------------------------------------------
# 2️⃣ Perform dynamic subsampling
# -------------------------------------------------------------
rule subsample_reads:
    input:
        os.path.join(config["tmp_dir"], "merged", "all_samples_trimmed.fastq")
    output:
        os.path.join(
            config["output_dir"],
            "03-subsample",
            "sample_size_{size}",
            "all_samples_subsampled_{size}.fastq"
        )
    log:
        os.path.join(config["log_dir"], "03-subsample", "{size}.log")
    params:
        sizes=get_dynamic_subsample_sizes
    container:
        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
    conda:
        "../envs/snakemake_usearch.yml"
    threads: 2
    resources:
        mem_mb=64000,     # 👈 give it 64 GB
        runtime=120,      # minutes or use SLURM default if not handled
    message:
        "Subsampling merged reads to size {wildcards.size}"
    shell:
        r"""
        exec &> "{log}"
        set -euxo pipefail

        mkdir -p "$(dirname {output})"
        input_file={input}
        size={wildcards.size}

        if [ "$size" = "all_reads" ]; then
            echo "Copying full dataset as all_reads subset"
            cp "$input_file" "{output}"
            exit 0
        fi

        nreads=$(grep -c '^+$' "$input_file" || true)
        if [ "$nreads" -lt "$size" ]; then
            echo "⚠️ Not enough reads ($nreads) for size $size, skipping..."
            : > "{output}"
            exit 0
        fi

        usearch -fastx_subsample "$input_file" \
            -sample_size $size \
            -fastqout {output} \
            -randseed 42 \
            -quiet
        """


rule merge_subsample_summaries:
    input:
        merged_fastq = os.path.join(config["tmp_dir"], "merged", "all_samples_trimmed.fastq"),
        subsamples = lambda wildcards: expand(
            os.path.join(
                config["output_dir"],
                "03-subsample",
                "sample_size_{size}",
                "all_samples_subsampled_{size}.fastq"
            ),
            size=get_dynamic_subsample_sizes(wildcards)
        )
    output:
        os.path.join(config["output_dir"], "03-subsample", "subsample_summary.tsv")
    log:
        os.path.join(config["log_dir"], "03-subsample", "merge_summary.log")
    params:
        replicate_list=",".join(sample_dirs)
    message:
        "Generating unified subsample summary file"
    shell:
        r"""
        exec &> "{log}"
        set -euxo pipefail

        merged_sample="all_samples"
        replicate_list="{params.replicate_list}"
        total_reads_file="$(realpath {input.merged_fastq})"

        if [ -s "$total_reads_file" ]; then
            total_reads=$(grep -c '^+$' "$total_reads_file" || echo 0)
        else
            echo "⚠️ Missing merged file: $total_reads_file"
            total_reads=0
        fi

        reads_after_trim=$total_reads  # identical here since merging uses trimmed reads

        echo -e "Merged_Sample\tReplicates\tTotal_Reads\tReads_After_Trim\tSubsample_Size\tReads_Subsampled" > {output}

        for subsample_file in {input.subsamples}; do
            size=$(basename "$subsample_file" | sed 's/all_samples_subsampled_//' | sed 's/.fastq//')
            if [ -s "$subsample_file" ]; then
                reads_subsampled=$(grep -c '^+$' "$subsample_file" || echo 0)
            else
                reads_subsampled=0
            fi
            echo -e "${{merged_sample}}\t${{replicate_list}}\t${{total_reads}}\t${{reads_after_trim}}\t${{size}}\t${{reads_subsampled}}" >> {output}
        done

        echo "✅ Summary written to {output}"
        """




