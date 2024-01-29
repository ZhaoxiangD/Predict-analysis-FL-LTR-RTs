#!/bin/bash
set -e
set -u
set -o pipefail
sample_name=$(basename -s '_finder_result.txt' "$1")
echo "${sample_name}" | python3 ../script/finder-result-process.py > ${sample_name}_result.xls

