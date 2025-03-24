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

view_alignment() {
    echo "BEGIN SAM CONVERSION $BASENAME"
    echo "INPUT BAM: ${WORKING_DIRECTORY}/data/alignments/${BASENAME}"

    # Check input BAM exists
    INPUT_BAM="${WORKING_DIRECTORY}/data/alignments/${BASENAME}"
    if [ ! -f "$INPUT_BAM" ]; then
        echo "Error: Input BAM file not found"
        return 1
    fi

    # Define output SAM path
    OUTPUT_SAM="${WORKING_DIRECTORY}/data/alignments/${BASENAME}.sam"

    echo "Starting BAM to SAM conversion..."
    samtools view "$INPUT_BAM" > "$OUTPUT_SAM"

    # Verify conversion success
    if [ ! -f "$OUTPUT_SAM" ]; then
        echo "Error: Failed to create SAM file"
        return 1
    fi

    echo "Successfully created: $OUTPUT_SAM"
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
        view_alignment
    done < "$SAMPLE_LIST"
}

main() {
    SAMPLE_LIST="$1"
    REFERENCE_GENOME="$2"
    WORKING_DIRECTORY="$(pwd)"
    
    process_samples
}

main "$@"
