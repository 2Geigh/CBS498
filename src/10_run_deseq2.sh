#!/bin/bash

run_deseq2() {
    WORKING_DIRECTORY="$(pwd)"

    Rscript ./src/10_run_deseq2.r
}

main() {
    run_deseq2 "$@"
}

main "$@"
