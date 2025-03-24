#!/bin/bash

visualize_data() {
    WORKING_DIRECTORY="$(pwd)"

    Rscript ./src/visualization.r
}

main() {
    visualize_data "$@"
}

main "$@"
