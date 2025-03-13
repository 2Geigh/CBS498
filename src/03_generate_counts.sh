# 03_generate_counts.sh
#!/bin/bash

generate_counts() {
    WORKING_DIRECTORY="$(pwd)"
    BAM_FILES=()
    
    for bam in "${WORKING_DIRECTORY}/data/alignments/"*.bam; do
        BAM_FILES+=("$bam")
    done
    
    featureCounts -a Arabidopsis_thaliana.TAIR10.60.gff3 \
        -o "${WORKING_DIRECTORY}/data/quant/output_data.txt" \
        "${BAM_FILES[@]}"
    
    cut -f1,7- "${WORKING_DIRECTORY}/data/quant/output_data.txt" > \
        "${WORKING_DIRECTORY}/data/quant/clean_counts.txt"
}

generate_counts "$@"
