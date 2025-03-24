#!/bin/bash

run_deseq2() {
    WORKING_DIRECTORY="$(pwd)"

    Rscript ./src/deseq2.r
}

main() {
    run_deseq2 "$@"
}

main "$@"
