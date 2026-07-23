#!/bin/bash
# Usage: bash add_lanes_to_filename.sh Sample_*/Lane_*/Unaligned

mkdir -p input_data

find "$@" -name '*.fq.gz' | while read -r f; do
  lane=$(sed -E 's|.*/Lane_([0-9]+)_.*|\1|' <<< "$f")
  ln -s "$(readlink -f "$f")" "input_data/L${lane}_$(basename "$f")"
done