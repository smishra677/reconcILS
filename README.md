# reconcILS

## Reconcile gene tree with species tree.

## Arguments
--spTree : Species tree in string
--gTree  : gene tree in string
--output : name of outputfile with extension .csv


## Output
Process,   Replicate, Gene_tree,   Species_Tree,    Duplication,   NNI,   Loss
reconcILS,     0,     "(B,(C,A));", "(A,(B,C));"    ,0,             1,     0



## Example Usage

python ./src/reconcILS.py --spTree '(A,(B,C));' --gTree '((A,C),B);' --output 'result.csv'