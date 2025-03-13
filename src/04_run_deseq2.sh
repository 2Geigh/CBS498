# 04_run_deseq2.sh
#!/bin/bash

run_deseq2() {
    WORKING_DIRECTORY="$(pwd)"
    Rscript << EOF
    # Existing R code remains unchanged
    # R automatically handles quoted paths from bash
EOF
}

run_deseq2 "$@"
