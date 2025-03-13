#!/bin/bash
# Enable exit on error for the entire script
set -e

# Make all scripts executable
chmod +x *.sh

# Run pipeline steps in order
./src/01_check_tools.sh
./src/02_run_pipeline.sh "$@"
./src/03_generate_counts.sh
./src/04_run_deseq2.sh