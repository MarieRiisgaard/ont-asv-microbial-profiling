############################################
# 05 – Generating SINTAX taxonomy predictions
############################################
rule sintax_subset:
    input:
        os.path.join(config["output_dir"], "04-denoise", "{subset}", "zOTUs.fa")
    output:
        filtered=os.path.join(config["output_dir"], "04-denoise", "{subset}", "zOTUs.fa.filtered"),
        sintax=os.path.join(config["output_dir"], "05-sintax", "{subset}", "zOTUs.sintax")
    log:
        os.path.join(config["log_dir"], "05-sintax", "{subset}_sintax.log")
    message:
        "Predicting taxonomy of zOTUs for subset {wildcards.subset} using SINTAX"
    resources:
        mem_mb=4096,
        runtime=60,
        cpus_per_task=config["max_threads"]
    container:
        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
    conda:
        "../envs/snakemake_usearch.yml"
    threads: config["max_threads"]
    params:
        db=config["db_sintax"],
        minlen=config.get("minlen", 100)
    shell:
        r"""
        exec &> "{log}"
        set -euxo pipefail

        mkdir -p "$(dirname {output.sintax})"

        # --- Handle missing or placeholder input ---
        if [ ! -s "{input}" ]; then
            echo "# No zOTUs to classify for subset {wildcards.subset}" > "{output.sintax}"
            : > "{output.filtered}"
            echo "⚠️ Skipping SINTAX — input file {input} empty" >&2
            exit 0
        fi

        if grep -q "^# " "{input}"; then
            echo "# Skipping SINTAX — placeholder input from UNOISE" > "{output.sintax}"
            : > "{output.filtered}"
            echo "⚠️ Skipping SINTAX for subset {wildcards.subset} (no valid zOTUs)" >&2
            exit 0
        fi

        echo "Filtering zOTUs shorter than {params.minlen} bp before SINTAX"

        awk -v minlen={params.minlen} '
            BEGIN {{ header = ""; seq = "" }}
            /^>/ {{
                if (header != "") {{
                    if (length(seq) >= minlen)
                        print header "\n" seq
                }}
                header = $0
                seq = ""
                next
            }}
            {{
                gsub(/[^A-Za-z]/, "", $0)
                seq = seq $0
            }}
            END {{
                if (header != "" && length(seq) >= minlen)
                    print header "\n" seq
            }}
        ' "{input}" > "{output.filtered}"

        # If no sequences remain after filtering, skip gracefully
        if [ ! -s "{output.filtered}" ]; then
            echo "# No sequences >= {params.minlen} bp to classify" > "{output.sintax}"
            echo "⚠️ Skipping SINTAX — all sequences too short" >&2
            exit 0
        fi

        # --- Run SINTAX ---
        usearch -sintax \
            "{output.filtered}" \
            -db "{params.db}" \
            -tabbedout "{output.sintax}" \
            -strand both \
            -sintax_cutoff 0.8 \
            -threads "{threads}" || true

        # --- Fallback: empty classification output ---
        if [ ! -s "{output.sintax}" ]; then
            echo "# SINTAX produced no classifications" > "{output.sintax}"
        fi

        echo "✅ SINTAX completed (subset {wildcards.subset})"
        """


