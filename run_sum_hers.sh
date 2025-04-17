#!/bin/bash

tagging_file="results/saliva_tagging.tagging"
results_folder="results/outputs"

while IFS=$'\t' read -r file_name ascertainment prevalence; do
    [[ "$file_name" == "file_name" ]] && continue

    if [[ -f "GWAS/${file_name}" ]]; then
        ldak --sum-hers "${results_folder}/output_${file_name}" \
             --tagfile "$tagging_file" \
             --summary "GWAS/${file_name}" \
             --ascertainment "$ascertainment" \
             --prevalence "$prevalence" \
             --check-sums NO
    else
        echo "File GWAS/${file_name} not found. Skipping..."
    fi
done < "GWAS_data.txt"
