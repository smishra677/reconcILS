# Copy Number Hemiplasy Simulator
Simulate gene trees under a model that allows copy number hemiplasy.

# Arguments
* stree
    + type = str
    + help = Path to file with newick formatted species tree to use for simulation. Branch lengths should be in coalescent units.

* mu_par
    + type = float
    + help = loss rate.

* lambda_par
    + type = float
    + help = duplication rate.

* reps
    + type = int
    + help = Number of gene families to simulate.

* output
    + type = str
    + help = Folder for storing results.

# Outputs

* rep_i.tre = final simulated gene family tree for replicate i.
* rep_i_alltrees.tre = all full subtrees (prior to the placement of mutations) for replicate i.
* rep_i.log = A log file indicating the total numbers of duplications and losses for the replicate. Also indicates when all copies were lost, and has counters for different types of deep coalescence events.
* rep_i.tsv = A table indicating for each branch in the species tree, which subtrees were concordant.
* rep_sptree_i.tre = An annotated species tree for matching branches to those in the tsv file.

# Example Usage
python3 simulator_v1g.py --stree sp_tree.tre --mu_par 0 --lambda_par 0.3 --reps 10 --output example