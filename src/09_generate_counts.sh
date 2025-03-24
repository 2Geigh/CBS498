#!/bin/bash

generate_counts() {
    WORKING_DIRECTORY="$(pwd)"
    OUTPUT_DESTINATION="${WORKING_DIRECTORY}/data/quant"
    TIMEPOINT_1_DIRECTORY="${OUTPUT_DESTINATION}/1h"
    TIMEPOINT_2_DIRECTORY="${OUTPUT_DESTINATION}/8h"
    TIMEPOINT_1_OUTPUT_DATA_PATH="${TIMEPOINT_1_DIRECTORY}/output_data.txt"
    TIMEPOINT_1_CLEAN_COUNTS_PATH="${TIMEPOINT_1_DIRECTORY}/clean_counts.txt"
    TIMEPOINT_2_OUTPUT_DATA_PATH="${TIMEPOINT_2_DIRECTORY}/output_data.txt"
    TIMEPOINT_2_CLEAN_COUNTS_PATH="${TIMEPOINT_2_DIRECTORY}/clean_counts.txt"

    TIMEPOINT_1_BAM_FILES=()
    TIMEPOINT_2_BAM_FILES=()

    REFERENCE_GENOME="Arabidopsis_thaliana.TAIR10.55.gtf"

    
    # Ensure the output directory and files exist
    mkdir -p "${OUTPUT_DESTINATION}"
    mkdir -p "${TIMEPOINT_1_DIRECTORY}"
    mkdir -p "${TIMEPOINT_2_DIRECTORY}"
    touch "${TIMEPOINT_1_OUTPUT_DATA_PATH}" "${TIMEPOINT_1_CLEAN_COUNTS_PATH}"
    
    for bam in "${WORKING_DIRECTORY}/data/alignments/"*.bam; do
        # Get the third character of the filename (excluding path)
        third_char="${bam##*/}"
        third_char="${third_char:2:1}"
        
        # Add to appropriate timepoint array based on third character
        case "$third_char" in
            0|1) TIMEPOINT_1_BAM_FILES+=("$bam") ;;
            *)   TIMEPOINT_2_BAM_FILES+=("$bam") ;;
        esac
    done
    
    featureCounts -p -a "${REFERENCE_GENOME}" \
        -o "${TIMEPOINT_1_OUTPUT_DATA_PATH}" \
        "${TIMEPOINT_1_BAM_FILES[@]}"
    
    cut -f1,7- "${TIMEPOINT_1_OUTPUT_DATA_PATH}" > \
        "${TIMEPOINT_1_CLEAN_COUNTS_PATH}"

    featureCounts -p -a "${REFERENCE_GENOME}" \
    -o "${TIMEPOINT_2_OUTPUT_DATA_PATH}" \
    "${TIMEPOINT_2_BAM_FILES[@]}"

    cut -f1,7- "${TIMEPOINT_2_OUTPUT_DATA_PATH}" > \
        "${TIMEPOINT_2_CLEAN_COUNTS_PATH}"
}

generate_counts "$@"
