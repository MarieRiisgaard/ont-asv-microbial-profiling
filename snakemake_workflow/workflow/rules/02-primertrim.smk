###############################################
# 02-primertrim.smk — per-barcode trimming
###############################################

rule trim_primers:
    input:
        os.path.join(
            config["tmp_dir"],
            "01-sample_prep",
            "{sample}",
            "{sample}_filtered_renamed.fastq",
        ),
    output:
        fastq_trimmed=os.path.join(
            config["tmp_dir"],
            "02-trim_primers",
            "{sample}",
            "{sample}_trimmed.fastq",
        ),
        reads_trimmed=os.path.join(
            config["tmp_dir"],
            "totalreads_trimmed",
            "{sample}_totaltrimmedreads.csv",
        ),
    log:
        os.path.join(config["log_dir"], "02-trim_primers", "trim_primers_{sample}.log"),
    params:
        primers=config["primers"],
    message:
        "Trimming primers for {wildcards.sample}",
    resources:
        mem_mb=2048,
        runtime=120,
        cpus_per_task=config["max_threads"],
    container:
        "docker://ghcr.io/kasperskytte/snakemake_usearch:main",
    conda:
        "../envs/snakemake_usearch.yml",
    threads: config["max_threads"],
    shell:
        r"""
        exec &> "{log}"
        set -euxo pipefail

        echo "*** Primer trimming for {wildcards.sample}"

        # --- Primer trimming with cutadapt ---
        cutadapt -g {params.primers} \
            --revcomp \
            -o {output.fastq_trimmed}.tmp \
            --discard-untrimmed {input} \
            -j {threads}

        # --- Remove any empty sequences (length = 0) ---
        awk 'NR%4==1 {{h=$0}} NR%4==2 {{s=$0}} NR%4==3 {{p=$0}} NR%4==0 {{q=$0; if(length(s)>0) print h"\n"s"\n"p"\n"q}}' \
            {output.fastq_trimmed}.tmp > {output.fastq_trimmed}

        rm {output.fastq_trimmed}.tmp

        # --- Check output file not empty ---
        #if [ ! -s "{output.fastq_trimmed}" ]; then
        #    echo "Output file {output.fastq_trimmed} is empty after filtering! Exiting."
        #    exit 1
        #fi

        # --- Check output file not empty ---
        if [ ! -s "{output.fastq_trimmed}" ]; then
            echo "No reads passed filtering. Creating empty FASTQ."
            : > "{output.fastq_trimmed}"
            num_reads=0
        else
            # Count reads after trimming & filtering
            num_reads=$(grep -c '^+$' "{output.fastq_trimmed}")
        fi

        echo "{wildcards.sample},$num_reads" > "{output.reads_trimmed}"

        echo "*** Primer trimming and cleanup completed for {wildcards.sample}"



        # --- Count reads after trimming & filtering ---
        #num_reads=$(grep -c '^+$' {output.fastq_trimmed})
        #echo "{wildcards.sample},$num_reads" > "{output.reads_trimmed}"

        #echo "*** Primer trimming and cleanup completed for {wildcards.sample}"
        #"""

###############################################
# End of file
###############################################
