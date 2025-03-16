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

#!/bin/bash

alignment() {
    echo "BEGIN HISAT2 $BASENAME ALIGNMENT"
    echo "REFERENCE GENOME: $REFERENCE_GENOME"

    # Check reference genome
    if [ ! -f "$REFERENCE_GENOME" ]; then
        echo "Error: Reference genome not found"
        return 1
    fi

    INDEX_DIR="indexes"
    mkdir -p "$INDEX_DIR"
    
    INDEX_BASENAME=$(basename "${REFERENCE_GENOME%.fa}")
    INDEX_BASENAME=${INDEX_BASENAME// /_}
    HISAT2_INDEX="$INDEX_DIR/$INDEX_BASENAME"

    hisat2-build "$REFERENCE_GENOME" "$HISAT2_INDEX"

    if [ ! -f "${HISAT2_INDEX}.1.ht2" ]; then
        echo "Error: Failed to build index"
        return 1
    fi

    echo "Starting alignment..."
    hisat2 -q -x "$HISAT2_INDEX" \
        -1 "$(printf '%q' "${WORKING_DIRECTORY}/data/trimmed_reads/${BASENAME}_R1_001.fastq.paired.R1.fastq")" \
        -2 "$(printf '%q' "${WORKING_DIRECTORY}/data/trimmed_reads/${BASENAME}_R1_001.fastq.paired.R2.fastq")" | \
        samtools sort -o "$(printf '%q' "${WORKING_DIRECTORY}/data/alignments/${BASENAME}.bam")"
}

process_samples() {
    echo "########################"
    echo "BEGIN PROCESSING SAMPLES"

    while IFS= read -r line; do
        # Add path escaping for sample files
        IFS=',' read -r READ1 READ2 <<< "$line"
        READ1="$READ1"
        READ2="$READ2"
        
        # Extract basename safely
        BASENAME=$(basename "$READ1" | sed 's/[._]R1.*//')
        BASENAME="${BASENAME// /_}"  # Replace spaces with underscores

        echo ""
        echo "############"
        echo "############"
        echo "Processing sample: $BASENAME"
        alignment
    done < "$SAMPLE_LIST"
}

main() {
    SAMPLE_LIST="$1"
    REFERENCE_GENOME="$2"
    WORKING_DIRECTORY="$(pwd)"
    
    process_samples
}

main "$@"
