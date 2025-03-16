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

QC_1() {
    echo "BEGIN FASTQC-ING $BASENAME"
        # QC with quoted outputs
        fastqc "$READ1" "$READ2" -o "${WORKING_DIRECTORY}/data/fastqc1/"
        echo "FASTQC-ING $BASENAME COMPLETE"
}

process_samples() {
    echo "########################"
    echo "BEGIN PROCESSING SAMPLES"
    while IFS= read -r line; do
        # Add path escaping for sample files
        IFS=',' read -r READ1 READ2 <<< "$line"
        
        # Extract basename safely
        BASENAME=$(basename "$READ1" | sed 's/[._]R1.*//')
        BASENAME="${BASENAME// /_}"  # Replace spaces with underscores

        echo " "
        echo " "
        echo " "
        echo " "
        echo " "
        echo "Processing sample: $BASENAME"

        QC_1
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
