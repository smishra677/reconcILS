#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

simphy="/N/slate/samishr/reconcILS_03_18/simphy/simphy"
base="/N/slate/samishr/reconcILS_03_18/EXPERIMENT_6/reconcILS/experiments"
#species_file="$base/species_40_scaled_10_simphy.tre"
species_file="$base/species_40_with_branch_lengths_10.tre"

i=0

while IFS= read -r species_tree || [[ -n "$species_tree" ]]; do
    species_tree="${species_tree%$'\r'}"
    [[ -z "${species_tree//[[:space:]]/}" ]] && continue

    i=$((i + 1))
    printf -v tree_id "%02d" "$i"
    output_dir="$base/simphy_40_base"

    "$simphy" \
        -S "$species_tree" \
        -LB f:0.00566 \
        -LD f:0.0045667 \
        -RL f:1000 \
        -SP f:10 \
        -OD 1 \
        -O "$output_dir" \
        -SG f:0.5

    cat "$output_dir"/1/g_trees*.trees > "$output_dir/g_trees.tre"

    echo "Completed: $output_dir"
done < "$species_file"