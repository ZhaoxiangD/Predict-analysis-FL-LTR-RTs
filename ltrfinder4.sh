#!/bin/bash
set -e
set -u
set -o pipefail
sample_name=$(basename -s '.fna' "$1")
ltr_finder "$1" -D 44000 -d 2500 -L 6000 -l 80 -p 20 -j 0.8 -S 6 -w2 > ${sample_name}_finder_result.txt

