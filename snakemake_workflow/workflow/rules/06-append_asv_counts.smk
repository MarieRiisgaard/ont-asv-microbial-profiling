############################################
# 06 – Append ASV counts and reference matches (exact + weak)
############################################
rule append_asv_counts:
    input:
        summary=os.path.join(config["output_dir"], "03-subsample", "subsample_summary.tsv"),
        zotus=expand(
            os.path.join(config["output_dir"], "04-denoise", "{subset}", "zOTUs.fa.filtered"),
            subset=subset_dirs
        )
    output:
        summary_with_asvs=os.path.join(config["output_dir"], "subsample_summary_with_ASVs.tsv"),
        detailed_matches=os.path.join(config["output_dir"], "subsample_summary_ASV_matches.tsv")
    log:
        os.path.join(config["log_dir"], "06-append_asv_counts.log")
    params:
        reference_amplicons=config.get("reference_amplicons", "")
    message:
        "Appending ASV counts and generating reference match report (exact + weak matches)"
    container:
        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
    conda:
        "../envs/snakemake_usearch.yml"
    shell:
        r"""
        exec &> "{log}"
        set -euxo pipefail

        mkdir -p "$(dirname {output.summary_with_asvs})"

        # --- If no summary file, stop gracefully ---
        if [ ! -s "{input.summary}" ]; then
            echo -e "Subset\tASVs\tASV_ref_matches\tAll_reads_ASVs\tAll_reads_ref_matches" > "{output.summary_with_asvs}"
            echo -e "Subset\tASV_ID\tReference_ID\tIdentity\tMatch_Type" > "{output.detailed_matches}"
            exit 0
        fi

        tmp_counts=$(mktemp)
        tmp_ref=$(mktemp)
        echo -e "Subset\tASVs" > "$tmp_counts"
        echo -e "Subset\tASV_ref_matches" > "$tmp_ref"
        echo -e "Subset\tASV_ID\tReference_ID\tIdentity\tMatch_Type" > "{output.detailed_matches}"

        # --- Collect ASV counts + perform matching ---
        if [ -n "{params.reference_amplicons}" ] && [ -s "{params.reference_amplicons}" ]; then
            echo "Reference file found: {params.reference_amplicons}"

            for zotus_file in {input.zotus}; do
                subset=$(basename "$(dirname "$zotus_file")")

                # Skip placeholder files
                if grep -q "^# " "$zotus_file"; then
                    echo -e "${{subset}}\t0" >> "$tmp_counts"
                    echo -e "${{subset}}\t0" >> "$tmp_ref"
                    continue
                fi

                # --- Count ASVs ---
                n_asv=$(grep -c '^>' "$zotus_file" || echo 0)
                echo -e "${{subset}}\t${{n_asv}}" >> "$tmp_counts"

                # --- Exact matches (1.00 identity) ---
                usearch -search_exact "$zotus_file" \
                    -db "{params.reference_amplicons}" \
                    -strand both \
                    -userout "${{zotus_file}}_exact.txt" \
                    -userfields query+target+id \
                    -quiet || true

                n_exact=$(grep -c . "${{zotus_file}}_exact.txt" || echo 0)
                echo -e "${{subset}}\t${{n_exact}}" >> "$tmp_ref"

                if [ -s "${{zotus_file}}_exact.txt" ]; then
                    awk -v subset="${{subset}}" '{{print subset"\t"$1"\t"$2"\t"$3"\tExact"}}' "${{zotus_file}}_exact.txt" >> "{output.detailed_matches}"
                fi

                # --- Weak / approximate matches (≥0.75 identity) ---
                usearch -usearch_global "$zotus_file" \
                    -db "{params.reference_amplicons}" \
                    -id 0.75 \
                    -strand both \
                    -userout "${{zotus_file}}_weak.txt" \
                    -userfields query+target+id \
                    -quiet \
                    -maxaccepts 0 -maxrejects 0 || true

                if [ -s "${{zotus_file}}_weak.txt" ]; then
                    awk -v subset="${{subset}}" '{{print subset"\t"$1"\t"$2"\t"$3"\tWeak"}}' "${{zotus_file}}_weak.txt" >> "{output.detailed_matches}"
                fi

                # --- Identify ASVs with no match (id < 0.75) ---
                all_asvs=$(grep '^>' "$zotus_file" | sed 's/^>//' || true)
                matched_asvs=$( (awk '{{print $1}}' "${{zotus_file}}_exact.txt" 2>/dev/null; awk '{{print $1}}' "${{zotus_file}}_weak.txt" 2>/dev/null) | sort -u )
                unmatched=$(comm -23 <(echo "$all_asvs" | sort) <(echo "$matched_asvs" | sort) || true)

                if [ -n "$unmatched" ]; then
                    while read -r asv; do
                        echo -e "${{subset}}\t${{asv}}\tNA\t0.00\tNo_match" >> "{output.detailed_matches}"
                    done <<< "$unmatched"
                fi
            done
        else
            echo "⚠️ No reference_amplicons specified or file missing — skipping identity checks."
            for zotus_file in {{input.zotus}}; do
                subset=$(basename "$(dirname "$zotus_file")")
                n_asv=$(grep -c '^>' "$zotus_file" || echo 0)
                echo -e "${{subset}}\t${{n_asv}}" >> "$tmp_counts"
                echo -e "${{subset}}\tNA" >> "$tmp_ref"
            done
        fi

        # --- Merge ASV counts + reference matches + all_reads data ---
        awk -v FS="\t" -v OFS="\t" -v cnt="$tmp_counts" -v ref="$tmp_ref" '
        BEGIN {{
            while ((getline line < cnt) > 0) {{
                if (line ~ /^Subset/) continue
                split(line, f, "\t")
                a[f[1]] = f[2]
            }}
            close(cnt)
            while ((getline line < ref) > 0) {{
                if (line ~ /^Subset/) continue
                split(line, f, "\t")
                r[f[1]] = f[2]
            }}
            close(ref)
            all_asv = (("all_reads" in a) ? a["all_reads"] : 0)
            all_ref = (("all_reads" in r) ? r["all_reads"] : "NA")
        }}
        NR == 1 {{
            print $0, "ASVs", "ASV_ref_matches", "All_reads_ASVs", "All_reads_ref_matches"
            next
        }}
        {{
            subset = $5  # column 5 = Subsample_Size
            asv = (subset in a ? a[subset] : 0)
            refv = (subset in r ? r[subset] : "NA")
            print $0, asv, refv, all_asv, all_ref
        }}' "{input.summary}" > "{output.summary_with_asvs}"

        rm "$tmp_counts" "$tmp_ref"
        echo "✅ append_asv_counts completed successfully."
        """



####    ########################################
# 06     – Append ASV counts and reference matches
####    ########################################
#rul    e append_asv_counts:
#        input:
#            summary=os.path.join(config["output_dir"], "03-subsample", "subsample_summary.tsv"),
#        zotus=expand(
#            os.path.join(config["output_dir"], "04-denoise", "{subset}", "zOTUs.fa.filtered"),
#            subset=subset_dirs
#        )
#    output:
#        summary_with_asvs=os.path.join(config["output_dir"], "subsample_summary_with_ASVs.tsv"),
#        detailed_matches=os.path.join(config["output_dir"], "subsample_summary_ASV_matches.tsv")
#    log:
#        os.path.join(config["log_dir"], "06-append_asv_counts.log")
#    params:
#        reference_amplicons=config.get("reference_amplicons", "")
#    message:
#        "Appending ASV counts and generating reference match report"
#    container:
#        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
#    conda:
#        "../envs/snakemake_usearch.yml"
#    shell:
#        r"""
#        exec &> "{log}"
#        set -euxo pipefail
#
#        mkdir -p $(dirname {output.summary_with_asvs})
#
#        # --- If no summary file, stop gracefully ---
#        if [ ! -s "{input.summary}" ]; then
#            echo -e "Subset\tASVs\tASV_ref_matches\tAll_reads_ASVs\tAll_reads_ref_matches" > "{output.summary_with_asvs}"
#            echo -e "Subset\tASV_ID\tReference_ID\tIdentity" > "{output.detailed_matches}"
#            exit 0
#        fi
#
#        # --- Collect ASV counts per subset ---
#        tmp_counts=$(mktemp)
#        echo -e "Subset\tASVs" > "$tmp_counts"
#        for zotus_file in {input.zotus}; do
#            subset=$(basename "$(dirname "$zotus_file")")
#            if [ -s "$zotus_file" ]; then
#                n_asv=$(grep -c '^>' "$zotus_file" || echo 0)
#            else
#                n_asv=0
#            fi
#            echo -e "${{subset}}\t${{n_asv}}" >> "$tmp_counts"
#        done
#
#        # --- Prepare reference match report ---
#        tmp_ref=$(mktemp)
#        echo -e "Subset\tASV_ref_matches" > "$tmp_ref"
#        echo -e "Subset\tASV_ID\tReference_ID\tIdentity" > "{output.detailed_matches}"
#
#        if [ -n "{params.reference_amplicons}" ] && [ -s "{params.reference_amplicons}" ]; then
#            echo "Reference file found: {params.reference_amplicons}"
#            for zotus_file in {input.zotus}; do
#                subset=$(basename "$(dirname "$zotus_file")")
#                if [ -s "$zotus_file" ]; then
#                    usearch -usearch_global "$zotus_file" \
#                        -db "{params.reference_amplicons}" \
#                        -id 0.95 \
#                        -strand both \
#                        -userout "${{zotus_file}}_refmatches.txt" \
#                        -userfields query+target+id \
#                        -quiet \
#                        -maxaccepts 0 -maxrejects 0 
#    
#                    n_match=$(awk '$3==1.0' "${{zotus_file}}_refmatches.txt" | wc -l || echo 0)
#
#                    if [ -s "${{zotus_file}}_refmatches.txt" ]; then
#                        awk -v subset="$subset" '{{{{print subset"\t"$1"\t"$2"\t"$3}}}}' "${{zotus_file}}_refmatches.txt" >> "{output.detailed_matches}"
#                    fi
#                else
#                    n_match=0
#                fi
#                echo -e "${{subset}}\t${{n_match}}" >> "$tmp_ref"
#            done
#        else
#            echo "⚠️ No reference_amplicons specified or file missing — skipping identity checks."
#            for zotus_file in {input.zotus}; do
#                subset=$(basename "$(dirname "$zotus_file")")
#                echo -e "${{subset}}\tNA" >> "$tmp_ref"
#            done
#        fi
#
#        # --- Merge ASV counts + reference matches + all_reads data ---
#        awk -v FS="\t" -v OFS="\t" -v cnt="$tmp_counts" -v ref="$tmp_ref" '
#        BEGIN {{{{
#            while ((getline line < cnt) > 0) {{{{
#                if (line ~ /^Subset/) continue
#                split(line, f, "\t")
#                a[f[1]] = f[2]
#            }}}}
#            close(cnt)
#            while ((getline line < ref) > 0) {{{{
#                if (line ~ /^Subset/) continue
#                split(line, f, "\t")
#                r[f[1]] = f[2]
#            }}}}
#            close(ref)
#
#            all_asv = (("all_reads" in a) ? a["all_reads"] : 0)
#            all_ref = (("all_reads" in r) ? r["all_reads"] : "NA")
#        }}}}
#        NR == 1 {{{{
#            print $0, "ASVs", "ASV_ref_matches", "All_reads_ASVs", "All_reads_ref_matches"
#            next
#        }}}}
#        /^barcode/ {{{{
#            subset = "sample_size_" $2
#            asv = (subset in a ? a[subset] : 0)
#            refv = (subset in r ? r[subset] : (("all_reads" in r) ? r["all_reads"] : "NA"))
#            print $0, asv, refv, all_asv, all_ref
#            next
#        }}}}
#        {{{{print}}}}
#        ' "{input.summary}" > "{output.summary_with_asvs}"
#
#        rm "$tmp_counts" "$tmp_ref"
#        echo "✅ append_asv_counts completed successfully."
#        """
#