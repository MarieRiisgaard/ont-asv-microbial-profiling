rule sample_prep:
    input:
        # function to list all fastq files per wildcard (subfolder/sample)
        # see https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#input-functions
        lambda wildcards: (
            glob.glob(
                os.path.join(config["input_dir"], wildcards.sample, "**", "*.fastq"),
                recursive=True,
            )
            + glob.glob(
                os.path.join(config["input_dir"], wildcards.sample, "**", "*.fq"),
                recursive=True,
            )
            + glob.glob(
                os.path.join(
                    config["input_dir"], wildcards.sample, "**", "*.fastq.gz"
                ),
                recursive=True,
            )
            + glob.glob(
                os.path.join(config["input_dir"], wildcards.sample, "**", "*.fq.gz"),
                recursive=True,
            )
        ),
    output:
        fastq=temp(
            os.path.join(
                config["tmp_dir"], "01-sample_prep", "{sample}", "{sample}.fastq"
            )
        ),
        total_reads_file=temp(
            os.path.join(config["tmp_dir"], "totalreads", "{sample}_totalreads.csv")
        ),
        sample_renamed=temp(
            os.path.join(
                config["tmp_dir"],
                "01-sample_prep",
                "{sample}",
                "{sample}_renamed.fastq",
            )
        ),
        fastq_filtered=temp(
            os.path.join(
                os.path.join(
                    config["tmp_dir"],
                    "01-sample_prep",
                    "{sample}",
                    "{sample}_filtered_renamed.fastq",
                )
            )
        ),
        totalfilteredreads_file=temp(
            os.path.join(
                config["tmp_dir"],
                "totalreads_filtered",
                "{sample}_totalfilteredreads.csv",
            )
        ),
    log:
        os.path.join(config["log_dir"], "01-sample_prep", "sample_prep_{sample}.log"),
    resources:
        mem_mb=lambda wc, input: max(5 * input.size_mb, 1024),
        runtime=30,
        cpus_per_task=1,
    container:
        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
    conda:
        "../envs/snakemake_usearch.yml"
    params:
        sample_sep=config["sample_sep"],
        filtlong_args=config["filtlong_args"],
    threads: 1
    message:
        "{wildcards.sample}: Filtering and preparing reads"
    shell:  
        """
        exec &> "{log}"
        set -euxo pipefail

        # Decompress only if compressed (.gz), otherwise just copy (EDITED FROM KSA version)
        echo "*** Decompressing and concatenating fastq files"
        if [[ "{input}" == *.gz ]]; then
            gunzip -cdfq {input} > {output.fastq}
        else
            cat {input} > {output.fastq}
        fi

        # calc total number of reads, TODO: count "+" lines instead
        echo "*** Calculating total number of reads before any filtering"
        num_reads=$(grep -c '^+$' {output.fastq})
        echo "{wildcards.sample},$num_reads" > "{output.total_reads_file}"

        echo "*** Renaming reads with sample name"
        usearch -fastx_relabel \
        "{output.fastq}" \
        -prefix "{wildcards.sample}{params.sample_sep}" \
        -fastqout "{output.sample_renamed}"

        echo "*** Filtering reads with Filtlong"
        filtlong {params.filtlong_args} {output.sample_renamed} > {output.fastq_filtered}

        # calc total number of reads
        echo "*** Calculating total number of reads after filtering"
        num_reads=$(grep -c '^+$' {output.fastq_filtered})
        echo "{wildcards.sample},$num_reads" > "{output.totalfilteredreads_file}"
        """


rule concatenate_total_reads_files:
    input:
        total_reads_file=expand(
            os.path.join(config["tmp_dir"], "totalreads", "{sample}_totalreads.csv"),
            sample=sample_dirs,
        ),
        totalfilteredreads_file=expand(
            os.path.join(
                config["tmp_dir"],
                "totalreads_filtered",
                "{sample}_totalfilteredreads.csv",
            ),
            sample=sample_dirs,
        ),
    output:
        total_reads_file=os.path.join(config["output_dir"], "totalreads.csv"),
        totalfilteredreads_file=os.path.join(
            config["output_dir"], "totalreads_filtered.csv"
        ),
    log:
        os.path.join(
            config["log_dir"], "01-sample_prep", "concatenate_total_reads_files.log"
        ),
    message:
        "Concatenating total reads files"
    resources:
        mem_mb=512,
        runtime=30,
        cpus_per_task=1,
    container:
        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
    conda:
        "../envs/snakemake_usearch.yml"
    threads: 1
    shell:
        """
        exec &> "{log}"
        set -euxo pipefail
        
        cat {input.total_reads_file} > {output.total_reads_file}
        cat {input.totalfilteredreads_file} > {output.totalfilteredreads_file}
        """
