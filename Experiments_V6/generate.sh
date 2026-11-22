# Experiment generation and analysis commands (V6).
# This is a record of individual commands, not a script to run end-to-end.
# Some commands use historical cluster paths or contain scratch/legacy syntax.
# Commands and parameters below are retained exactly as originally written.
# generate.r writes ten new 40-tip trees and a 10x SimPhy-scaled counterpart.
# The last SimPhy loop reads the scaled tree file (it does not scale trees itself).

# Simulate 1,000 gene trees on the 40-species tree at the baseline duplication/loss rates.
python3  /N/slate/samishr/reconcILS_03_18/dupcoal_12/dupcoal-main/dupcoal_v3.0.0.py --stree species_40_with_branch_lengths.tre --mu_par 0.03 --lambda_par 0.03 --reps 1000 --output Data_with_duplication_40_species_0.03_0.03
# Repeat with a lower duplication rate (lambda = 0.01).
python3 /N/slate/samishr/reconcILS_03_18/dupcoal_12/dupcoal-main/dupcoal_v3.0.0.py --stree species_40_with_branch_lengths.tre --mu_par 0.03 --lambda_par 0.01 --reps 1000 --output Data_with_duplication_40_species_0.03_0.01
# Repeat with a lower loss rate (mu = 0.01).
python3 /N/slate/samishr/reconcILS_03_18/dupcoal_12/dupcoal-main/dupcoal_v3.0.0.py --stree species_40_with_branch_lengths.tre --mu_par 0.01 --lambda_par 0.03 --reps 1000 --output Data_with_duplication_40_species_0.01_0.03

# Generate a three-species control with no duplications or losses (ILS only).
python3 /N/slate/samishr/reconcILS_03_18/dupcoal_12/dupcoal-main/dupcoal_v3.0.0.py --stree 3sp.tre --mu_par 0 --lambda_par 0 --reps 1000 --output Data_with_duplication_3_species_0_0


# Reconcile the baseline 40-species dupcoal dataset with default event costs.
python test_reconcILS_default_40.py \
  --data Data_with_duplication_40_species_0.03_0.03 \
  --species-tree species_40_with_branch_lengths.tre \
  --duplication-cost 1.1 \
  --ils-cost 1 \
  --loss-cost 1 \
  --workers 20

# Scale branch lengths by 100 using scale.py.
# This produces a _100.tre file; the later x1000 simulation requires a separate _1000.tre file.
python ./scale.py --input species_40_with_branch_lengths.tre --scale 100


# Simulate on the long-branch tree with correspondingly reduced event rates.
python3  /N/slate/samishr/reconcILS_03_18/dupcoal_12/dupcoal-main/dupcoal_v3.0.0.py --stree species_40_with_branch_lengths_1000.tre --mu_par 0.00003 --lambda_par 0.00003 --reps 1000 --output Data_with_duplication_40_species_0.00003_0.00003_x1000

# Generate fungal gene trees using the scaled fungal species tree.
python3 /N/slate/samishr/reconcILS_03_18/dupcoal_12/dupcoal-main/dupcoal_v3.0.0.py --stree fungi_scaled_0.1.tre --mu_par 0.0015 --lambda_par 0.003 --reps 1000 --output Data_with_duplication__fungi_1_0_1_0.0015_0.003


# Generate a three-species dataset with higher duplication and loss rates.
python3 /N/slate/samishr/reconcILS_03_18/dupcoal_12/dupcoal-main/dupcoal_v3.0.0.py --stree 3sp.tre --mu_par 0.15 --lambda_par 0.3 --reps 1000 --output Data_with_duplication_3_0.15_0.3

# Dataset name noted for reference (this line is not a runnable command).
Data_with_duplication_3_0.15_0.3


# Repeat the fungal simulation command above.
python3 /N/slate/samishr/reconcILS_03_18/dupcoal_12/dupcoal-main/dupcoal_v3.0.0.py --stree fungi_scaled_0.1.tre --mu_par 0.0015 --lambda_par 0.003 --reps 1000 --output Data_with_duplication__fungi_1_0_1_0.0015_0.003
# Repeat the three-species simulation at mu = 0.15, lambda = 0.3.
python3 /N/slate/samishr/reconcILS_03_18/dupcoal_12/dupcoal-main/dupcoal_v3.0.0.py --stree 3sp.tre --mu_par 0.15 --lambda_par 0.3 --reps 1000 --output Data_with_duplication_3_0.15_0.3

# Generate a second three-species dataset with mu = 0.1, lambda = 0.3.
python3 /N/slate/samishr/reconcILS_03_18/dupcoal_12/dupcoal-main/dupcoal_v3.0.0.py --stree 3sp.tre --mu_par 0.1 --lambda_par 0.3 --reps 1000 --output Data_with_duplication_3_0.1_0.3



# Compare reconcILS, ETE3, and DLCpar on the first three-species dataset.
python3 trigger_dlcpar.py   --species-tree 3sp.tre   --duplication-cost 1.1   --ils-cost 1   --loss-cost 1   --workers 20   --data Data_with_duplication_3_0.15_0.3  --dlcpar-python /N/slate/samishr/reconcILS_03_18/dlcpar-master/py2env/bin/python
# Run the same comparison for the second three-species dataset.
python3 trigger_dlcpar.py   --species-tree 3sp.tre   --duplication-cost 1.1   --ils-cost 1   --loss-cost 1   --workers 20   --data Data_with_duplication_3_0.1_0.3  --dlcpar-python /N/slate/samishr/reconcILS_03_18/dlcpar-master/py2env/bin/python


# Compare methods on fungal data using alternative event costs (2, 2, 1).
python3 trigger_dlcpar.py   --species-tree fungi_scaled_0.1.tre   --duplication-cost 2   --ils-cost 1   --loss-cost 2   --workers 20   --data Data_with_duplication__fungi_1_0_1_0.0015_0.003   --dlcpar-python /N/slate/samishr/reconcILS_03_18/dlcpar-master/py2env/bin/python
# Repeat the fungal comparison with default costs (1.1, 1, 1).
python3 trigger_dlcpar.py   --species-tree fungi_scaled_0.1.tre   --duplication-cost 1.1   --ils-cost 1   --loss-cost 1 --workers 20   --data Data_with_duplication__fungi_1_0_1_0.0015_0.003   --dlcpar-python /N/slate/samishr/reconcILS_03_18/dlcpar-master/py2env/bin/python


# Long-branch dataset name noted for reference (not a shell command).
Data_with_duplication_40_species_0.00003_0.00003_x1000

# Reconcile the long-branch 40-species dataset.
python test_reconcILS_default_40.py \
  --data Data_with_duplication_40_species_0.00003_0.00003_x1000 \
  --species-tree species_40_with_branch_lengths.tre \
  --duplication-cost 1.1 \
  --ils-cost 1 \
  --loss-cost 1 \
  --workers 20
# Historical alternative reconciliation commands below.
# NOTE: Python-style triple quotes are not Bash comments; this block is kept as originally written.
'''
# Long-branch comparison repeated with default costs.
python test_reconcILS_default_40.py \
  --data Data_with_duplication_40_species_0.00003_0.00003_x1000 \
  --species-tree species_40_with_branch_lengths.tre \
  --duplication-cost 1.1 \
  --ils-cost 1 \
  --loss-cost 1 \
  --workers 20

# Reconcile the three-species ILS-only control.
python test_reconcILS_default_40.py \
  --data Data_with_duplication_3_species_0_0 \
  --species-tree 3sp.tre \
  --duplication-cost 1.1 \
  --ils-cost 1 \
  --loss-cost 1 \
  --workers 20


python test_reconcILS_default_40.py \
# Baseline 40-species reconciliation at default costs.
  --data Data_with_duplication_40_species_0.03_0.03 \
  --species-tree species_40_with_branch_lengths.tre \
  --duplication-cost 1.1 \
  --ils-cost 1 \
  --loss-cost 1 \
  --workers 20

# Test a higher loss cost (duplication 3, loss 4, ILS 1).
python test_reconcILS_default_40.py \
  --data Data_with_duplication_40_species_0.03_0.03 \
  --species-tree species_40_with_branch_lengths.tre \
  --duplication-cost 3 \
  --ils-cost 1 \
  --loss-cost 4 \
  --workers 20

# Test a higher duplication cost (duplication 4, loss 3, ILS 1).
  python test_reconcILS_default_40.py \
  --data Data_with_duplication_40_species_0.03_0.03 \
  --species-tree species_40_with_branch_lengths.tre \
  --duplication-cost 4 \
  --ils-cost 1 \
  --loss-cost 3 \
  --workers 20

# Test a higher ILS cost (duplication 1, loss 1, ILS 5).
python test_reconcILS_default_40.py \
  --data Data_with_duplication_40_species_0.03_0.03 \
  --species-tree species_40_with_branch_lengths.tre \
  --duplication-cost 1 \
  --ils-cost 5 \
  --loss-cost 1 \
  --workers 20
'''

# Process ten numbered 40-species datasets with default reconciliation costs.
for i in {01..10}; do
    echo "$i"
    python test_reconcILS_default_40.py \
      --data "Data_with_duplication_40_species_0.03_0.03_$i" \
      --species-tree species_40_with_branch_lengths.tre \
      --duplication-cost 1.1 \
      --ils-cost 1 \
      --loss-cost 1 \
      --workers 20
done


# Reconcile the 40-species ILS-only dataset.
python test_reconcILS_default_40.py \
  --data   Data_with_duplication_40_species_just_ILS \
  --species-tree species_40_with_branch_lengths.tre \
  --duplication-cost 1.1 \
  --ils-cost 1 \
  --loss-cost 1 \
  --workers 20


# Reconcile a dataset with low loss rate and higher duplication rate.
python test_reconcILS_default_40.py \
  --data   Data_with_duplication_40_species_0.0001_0.03 \
  --species-tree species_40_with_branch_lengths.tre \
  --duplication-cost 1.1 \
  --ils-cost 1 \
  --loss-cost 1 \
  --workers 20

# Reconcile a dataset with higher loss rate and low duplication rate.
python test_reconcILS_default_40.py \
  --data   Data_with_duplication_40_species_0.03_0.0001 \
  --species-tree species_40_with_branch_lengths.tre \
  --duplication-cost 1.1 \
  --ils-cost 1 \
  --loss-cost 1 \
  --workers 20


# Run the SimPhy reconciliation trigger (legacy command; the slash sequence needs checking).
python trigger_simphy.py \ \
  --workers 20


# Generate the low-loss, high-duplication dupcoal dataset used above.
python3  /N/slate/samishr/reconcILS_03_18/dupcoal_12/dupcoal-main/dupcoal_v3.0.0.py --stree species_40_with_branch_lengths.tre --mu_par 0.0001 --lambda_par 0.03 --reps 1000 --output Data_with_duplication_40_species_0.0001_0.03


# Generate the high-loss, low-duplication dupcoal dataset used above.
python3  /N/slate/samishr/reconcILS_03_18/dupcoal_12/dupcoal-main/dupcoal_v3.0.0.py --stree species_40_with_branch_lengths.tre --mu_par 0.03 --lambda_par 0.0001 --reps 1000 --output Data_with_duplication_40_species_0.03_0.0001

# Separate shell-loop recipe begins here; kept in its original position.
#!/bin/bash
#!/bin/bash
set -euo pipefail

# Read one species tree at a time from the ten-tree file.
TREE_FILE="species_40_scaled_10.tre"
TEMP_TREE="species_tree_current.tre"

# Loop over ten independently scaled species-tree datasets.
for i in {01..10}; do
    echo "Running dataset $i"

    line_num=$((10#$i))

    # Extract the i-th species tree into a temporary file
    sed -n "${line_num}p" "$TREE_FILE" > "$TEMP_TREE"

    if [[ ! -s "$TEMP_TREE" ]]; then
        echo "Error: No species tree found on line $line_num"
        exit 1
    fi

# Reconcile each numbered dupcoal dataset against its corresponding species tree.
    python3 test_reconcILS_default_40.py \
        --data "Data_with_duplication_40_species_0.03_0.03_$i" \
        --species-tree "$TEMP_TREE" \
        --duplication-cost 1.1 \
        --ils-cost 1 \
        --loss-cost 1 \
        --workers 20

    echo "Completed dataset $i"
done

rm -f "$TEMP_TREE"
echo "All 10 datasets completed."
# Additional baseline reconciliation command (dataset name ends with an underscore).
python3 test_reconcILS_default_40.py \
    --data "Data_with_duplication_40_species_0.03_0.03_" \
    --species-tree species_40_with_branch_lengths.tre \
    --duplication-cost 1.1 \
    --ils-cost 1 \
    --loss-cost 1 \
    --workers 20


# Introduce gene-tree error with a mean of 5 NNI changes per tree.
python simulate_noise.py \
  --file /N/slate/samishr/reconcILS_03_18/EXPERIMENT_6/reconcILS/experiments/Data_with_duplication_40_species_0.03_0.03/ \
  --error 5 \


# Introduce gene-tree error with a mean of 10 NNI changes per tree.
python simulate_noise.py \
  --file /N/slate/samishr/reconcILS_03_18/EXPERIMENT_6/reconcILS/experiments/Data_with_duplication_40_species_0.03_0.03/ \
  --error 10 \


# Introduce gene-tree error with a mean of 15 NNI changes per tree.
python simulate_noise.py \
  --file /N/slate/samishr/reconcILS_03_18/EXPERIMENT_6/reconcILS/experiments/Data_with_duplication_40_species_0.03_0.03/ \
  --error 15 \



# Scratch Python argument-parser snippet; reference only, not Bash syntax.
  def parse1():
    parser = argparse.ArgumentParser(description="Add noise ")
    parser.add_argument('--file', type=str, help="Path to gene tree folder")
    parser.add_argument('--error', type=str, help="Error")
    args= parser.parse_args()



# Reconcile the trees with 5 NNI errors and save their results separately.
python3 test_reconcILS_default_40.py \
    --data "Data_with_duplication_40_species_0.03_0.03" \
    --species-tree species_40_with_branch_lengths.tre \
    --file-name "error_trees_5.tre" \
    --output-dir "Data_with_duplication_40_species_0.03_0.03\results_err_5" \
    --duplication-cost 1.1 \
    --ils-cost 1 \
    --loss-cost 1 \
    --workers 20

# Reconcile the trees with 10 NNI errors and save their results separately.
python3 test_reconcILS_default_40.py \
    --data "Data_with_duplication_40_species_0.03_0.03" \
    --species-tree species_40_with_branch_lengths.tre \
    --file-name "error_trees_10.tre" \
     --output-dir "Data_with_duplication_40_species_0.03_0.03\results_err_10" \
    --duplication-cost 1.1 \
    --ils-cost 1 \
    --loss-cost 1 \
    --workers 20


# All-paralogy is the default: python3 run_all_paralogy.py (no arguments).
# For SCO, provide the SCO gene tree instead.
# Run the single-copy ortholog (SCO) primate dataset.
python3 run_all_paralogy.py \
    --genes "./ZeroCol_ASTRAL_ML_SCO_MIN27.tre" \
     --output-dir "sco_efficient" \





# Separate SimPhy comparison loop begins here.
# Its species_40_scaled_10_simphy.tre input must be located in the working directory.
# The bundled copy is under "all_species tree/"; generate.r writes a new copy in "generated/".
#!/bin/bash
#!/bin/bash
set -euo pipefail

# Read the ten SimPhy-scaled species trees from a multi-tree file.
TREE_FILE="species_40_scaled_10_simphy.tre"
TEMP_TREE="species_tree_current.tre"

# Loop over the ten SimPhy datasets.
for i in {01..10}; do
    echo "Running dataset $i"

    line_num=$((10#$i))

    # Extract the i-th species tree into a temporary file
    sed -n "${line_num}p" "$TREE_FILE" > "$TEMP_TREE"

    if [[ ! -s "$TEMP_TREE" ]]; then
        echo "Error: No species tree found on line $line_num"
        exit 1
    fi

# Run reconcILS on each SimPhy dataset using its corresponding species tree.
    python3 trigger_simphy.py \
        --data "simphy_40_0045666666667_0045666666667_$i" \
        --species "$TEMP_TREE" \
        --out "simphy_40_0045666666667_0045666666667_${i}_results" \
        --workers 20

    echo "Completed dataset $i"
done

rm -f "$TEMP_TREE"
echo "All 10 datasets completed."


