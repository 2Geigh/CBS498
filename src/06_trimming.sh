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

trimming() {
    echo "BEGIN TRIMMOMATIC-ING $BASENAME"
    # Trimmomatic with full quoting
    #PATH_TO_TRIMMOMATIC="/home/twogeigh/miniconda3/share/trimmomatic-0.39-2/trimmomatic.jar" #Linux
    PATH_TO_TRIMMOMATIC="/Users/2geigh/miniforge3/share/trimmomatic-0.39-2/trimmomatic.jar" #MacOS

    java -jar "${PATH_TO_TRIMMOMATIC}" PE \
        "$READ1" "$READ2" \
        "${WORKING_DIRECTORY}/data/trimmed_reads/${BASENAME}.paired.R1.fastq" \
        "${WORKING_DIRECTORY}/data/trimmed_reads/${BASENAME}.unpaired.R1.fastq" \
        "${WORKING_DIRECTORY}/data/trimmed_reads/${BASENAME}.paired.R2.fastq" \
        "${WORKING_DIRECTORY}/data/trimmed_reads/${BASENAME}.unpaired.R2.fastq" \
        TRAILING:10 MINLEN:36
    echo "TRIMMOMATIC-ING $BASENAME COMPLETE"
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
        BASENAME=$(basename "$READ1" | sed 's/[._]R[[:digit:]]\+.*//')
        BASENAME="${BASENAME// /_}"  # Replace spaces with underscores

        echo " "
        echo " "
        echo " "
        echo " "
        echo " "
        echo "Processing sample: $BASENAME"
        
        trimming
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
