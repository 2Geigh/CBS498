# 03_generate_counts.sh
#!/bin/bash

generate_counts() {
    WORKING_DIRECTORY="$(pwd)"
    BAM_FILES=()
    
    # Ensure the directories exist
    mkdir -p "${WORKING_DIRECTORY}/data/quant"
    
    for bam in "${WORKING_DIRECTORY}/data/alignments/"*.bam; do
        BAM_FILES+=("$bam")
    done
    
    # Initialize output files if they do not exist
    touch "${WORKING_DIRECTORY}/data/quant/output_data.txt" "${WORKING_DIRECTORY}/data/quant/clean_counts.txt"
    
    featureCounts -p -a ARaport11_GTF_genes_transposons.current.gtf \
        -o "${WORKING_DIRECTORY}/data/quant/output_data.txt" \
        "${BAM_FILES[@]}"
    
    cut -f1,7- "${WORKING_DIRECTORY}/data/quant/output_data.txt" > \
        "${WORKING_DIRECTORY}/data/quant/clean_counts.txt"
}

generate_counts "$@"
