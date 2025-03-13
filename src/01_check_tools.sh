# 01_check_tools.sh
#!/bin/bash

check_tools() {
    echo "Checking required tools..."
    TOOLS=(fastqc trimmomatic hisat2 samtools featureCounts R)
    
    for tool in "${TOOLS[@]}"; do
        if ! command -v "$tool" &>/dev/null; then
            echo "Error: $tool is not installed"
            exit 1
        fi
    done
    echo "All tools are properly installed."
}

check_tools "$@"
