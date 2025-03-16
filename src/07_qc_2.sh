# 02_run_pipeline.sh
#!/bin/bash

WORKING_DIRECTORY="$(printf '%q' "$(pwd)")"

validate_arguments() {
    # @description Validate command line arguments
    # @returns 0 on success, exits script on failure

    echo "Number of arguments received: $#"
    echo "Arguments received: $@"

    if [ $# -ne 2 ]; then
        echo "Usage: <sample_list.txt> <reference_genome.fa>"
        exit 1
    fi
}

QC_2() {
    # Verify fastqc2 directory exists before running FastQC
    FASTQC2_DIR="${WORKING_DIRECTORY}/data/fastqc2"
    if [ ! -d "$FASTQC2_DIR" ]; then
        mkdir -p "$FASTQC2_DIR"
        if [ ! -d "$FASTQC2_DIR" ]; then
            echo "Error: Could not create FastQC output directory: $FASTQC2_DIR"
            exit 1
        fi
    fi
    
    echo "BEGIN FASTQC-ING TRIMMED $BASENAME"
    fastqc \
        "${WORKING_DIRECTORY}/data/trimmed_reads/${BASENAME}_R1_001.fastq.paired.R1.fastq" \
        "${WORKING_DIRECTORY}/data/trimmed_reads/${BASENAME}_R1_001.fastq.paired.R2.fastq" \
        -o "$FASTQC2_DIR"
    echo "FASTQC-ING TRIMMED $BASENAME COMPLETE"
}

process_samples() {
    echo "########################"
    echo "BEGIN PROCESSING SAMPLES"
    while IFS= read -r line; do
        # Add path escaping for sample files
        IFS=',' read -r READ1 READ2 <<< "$line"
        READ1="$(echo "$READ1" | sed 's/ /\\ /g')"
        READ2="$(echo "$READ2" | sed 's/ /\\ /g')"
        
        # Extract basename safely
        BASENAME=$(basename "$READ1" | sed 's/[._]R1.*//')
        BASENAME="${BASENAME// /_}"  # Replace spaces with underscores

        echo " "
        echo " "
        echo " "
        echo " "
        echo " "
        echo "Processing sample: $BASENAME"

        QC_2
    done < "$SAMPLE_LIST"
}


main() {
    @description Main program entry point
    validate_arguments "$@"
    SAMPLE_LIST="$1"
    REFERENCE_GENOME="$2"
    WORKING_DIRECTORY="$(pwd)"
    
    echo "Sample list: $SAMPLE_LIST"
    echo "Reference genome: $REFERENCE_GENOME"

    process_samples
}

main "$@"
