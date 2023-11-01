# reconcILS

## Reconcile gene tree with species tree.

# Installation 

## pip install .


## Arguments
--spTree : Species tree in string
--gTree  : gene tree in string
--output : name of outputfile with extension .csv
--D:  Duplication Cost
--L:  Loss Cost
--I:  ILS Cost
--V: Verbose Mode

## Output
Process,   Replicate,  Gene_tree,   Species_Tree,    Duplication,   NNI,   Loss
reconcILS,     0,     "(((B,C_1)C_1 1I+1D,C_2),A);", "(A,(B,C));"    ,0,             1,     0



## Example Usage

python ./src/reconcILS.py --spTree '(A,(B,C));' --gTree '((A,C),B);' --output 'result.csv'
