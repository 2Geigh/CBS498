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

verify_sample_list() {
    # @description Check if sample list file exists
    # @param $1 sample list path
    # @returns 0 on success, exits script on failure
    if [ ! -f "$1" ]; then
        echo "Error: Sample list file not found"
        exit 1
    fi
}

verify_reference_genome() {
    # @description Check if reference genome exists and is accessible
    # @param $1 reference genome path
    # @returns 0 on success, exits script on failure
    local reference_genome="$1"
    
    if [ -z "$reference_genome" ]; then
        echo "Error: Reference genome not provided."
        exit 1
    fi
    
    if [ ! -f "$reference_genome" ]; then
        echo "Error: Reference genome not found at '$reference_genome'"
        exit 1
    fi
}

setup_directories() {
    # @description Create and verify all required directories
    # @param $1 working directory path
    # @returns 0 on success, exits script on failure
    local working_dir="$1"
    
    local dirs=(
        "${working_dir}/data/fastqc1"
        "${working_dir}/data/trimmed_reads"
        "${working_dir}/data/alignments"
        "${working_dir}/data/quant"
        "${working_dir}/data/multiqc1"
        "${working_dir}/data/multiqc2"
        "${working_dir}/results"
    )
    
    echo "\nBEGIN DIRECTORY CHECK"
    
    for dir in "${dirs[@]}"; do
        mkdir -p "$dir"
        
        if [ ! -d "$dir" ]; then
            echo "Error: Could not create directory: $dir"
            exit 1
        fi
        
        if [ ! -w "$dir" ]; then
            echo "Error: No write permissions for directory: $dir"
            exit 1
        fi
    done
    
    echo "DIRECTORY CHECK COMPLETE"
}

QC_1() {
    echo "BEGIN FASTQC-ING $BASENAME"
        # QC with quoted outputs
        fastqc "$READ1" "$READ2" -o "${WORKING_DIRECTORY}/data/fastqc1/"
        echo "FASTQC-ING $BASENAME COMPLETE"
}

trimming() {
    echo "BEGIN TRIMMOMATIC-ING $BASENAME"
    # Trimmomatic with full quoting
    PATH_TO_TRIMMOMATIC="/home/twogeigh/miniconda3/share/trimmomatic-0.39-2/trimmomatic.jar"
    java -jar "${PATH_TO_TRIMMOMATIC}" PE \
        "$READ1" "$READ2" \
        "${WORKING_DIRECTORY}/data/trimmed_reads/${BASENAME}.paired.R1.fastq" \
        "${WORKING_DIRECTORY}/data/trimmed_reads/${BASENAME}.unpaired.R1.fastq" \
        "${WORKING_DIRECTORY}/data/trimmed_reads/${BASENAME}.paired.R2.fastq" \
        "${WORKING_DIRECTORY}/data/trimmed_reads/${BASENAME}.unpaired.R2.fastq" \
        TRAILING:10 MINLEN:36
    echo "TRIMMOMATIC-ING $BASENAME COMPLETE"
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
        "${WORKING_DIRECTORY}/data/trimmed_reads/${BASENAME}.paired.R1.fastq" \
        "${WORKING_DIRECTORY}/data/trimmed_reads/${BASENAME}.paired.R2.fastq" \
        -o "$FASTQC2_DIR"
    echo "FASTQC-ING TRIMMED $BASENAME COMPLETE"
}

alignment() {
    echo "BEGIN HISAT2 $BASENAME ALIGNMENT"
    echo "REFERENCE GENOME: $REFERENCE_GENOME"
    
    # Check reference genome
    [ -f "$REFERENCE_GENOME" ] || { echo "Error: Reference genome not found"; exit 1; }
    
    # Create index directory
    INDEX_DIR="indexes"
    mkdir -p "$INDEX_DIR"
    
    # Build index in the index directory
    INDEX_BASENAME=$(basename "${REFERENCE_GENOME%.fa}" | tr ' ' '_')
    HISAT2_INDEX="$INDEX_DIR/$INDEX_BASENAME"
    hisat2-build "$REFERENCE_GENOME" "$HISAT2_INDEX"
    
    # Verify index
    [ -f "${HISAT2_INDEX}.1.ht2" ] || { echo "Error: Failed to build index"; exit 1; }
    
    # Align & sort
    echo "Starting alignment..."
    hisat2 -q -x "$HISAT2_INDEX" \
           -1 "${WORKING_DIRECTORY}/data/trimmed_reads/${BASENAME}.paired.R1.fastq" \
           -2 "${WORKING_DIRECTORY}/data/trimmed_reads/${BASENAME}.paired.R2.fastq" | \
           samtools sort -o "${WORKING_DIRECTORY}/data/alignments/${BASENAME}.bam"
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
        BASENAME=$(basename "$READ1" | sed 's/[._]R1.*//')  # Adjust pattern for your filenames
        BASENAME="${BASENAME// /_}"  # Replace spaces with underscores

        echo " "
        echo " "
        echo " "
        echo " "
        echo " "
        echo "Processing sample: $BASENAME"

        QC_1
        
        trimming

        QC_2

        alignment
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
    
    verify_sample_list "$SAMPLE_LIST"

    verify_reference_genome "$REFERENCE_GENOME"
    
    setup_directories "$WORKING_DIRECTORY"

    process_samples
}

main "$@"
