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
* rep_i_alltrees_mutated.tre = all full subtrees (after the placement of mutations) for replicate i. Note that when trees only have a single branch, they do not show up appropriately in FigTree.
* rep_i_alltrees_unmutated.tre = all full subtrees (prior to the placement of mutations) for replicate i.
* rep_i.log = A log file indicating the total numbers of duplications and losses for the replicate. Also indicates when all copies were lost, and has counters for different types of deep coalescence events.
* rep_i.tsv = A table indicating for each branch in the species tree, which subtrees were concordant.
* rep_sptree_i.tre = An annotated species tree for matching branches to those in the tsv file.

# Log file explanation

* Total duplications: number of duplications
* Total losses: number of losses
* Deep Coalescence in the parent tree: number of deep coalescences in the parent tree
* Deep Coalescence in the daughter trees: number of deep coalescences in the daughter tree
* Copy Number Hemiplasy: number of cases in which a mutation is placed on a branch that does not exist in the branch of the species tree on which the mutation occurred.
* Deep Coalescence due to placing mutations: number of cases in which a mutation is placed on a branch that does exist in the species tree, but does not match the current branch of the species tree.
* Deep Coalescences when joining daughter trees: number of deep coalescences that occur when joining parent and daughter subtrees.
* Deep Coalescences when joining daughter trees (only count once for each subtree, DLCPar counting): number of deep coalescences that occur when joining parent and daughter subtrees. In this case, the maximum count for any subtree is 1. ILS either happened or did not. This is somewhat redundant with the above counter, except that it will be lower in some cases.
* Deep Coalescences due only to branch mismatch joining daughter trees: number of deep coalescenses that occur when joining happens in the first available branch of the species tree, but is to a branch of another gene tree that does not match that branch of the species tree in terms of taxon composition.

# Example Usage
python3 simulator_v1g.py --stree sp_tree.tre --mu_par 0 --lambda_par 0.3 --reps 10 --output example