# reconcILS

## Reconcile Gene Tree with Species Tree

## Requirements
This package has the following requirements:

- Python (3.x)
- Pandas
- Matplotlib
- uuid
- igraph




## Installation

You can install reconcILS using pip:

```bash
pip install .
```

## Arguments

- `--spTree`: Species tree in string.
- `--gTree`: Gene tree in string.
- `--output`: Name of the output file with the extension .csv.
- `--D`: Duplication Cost.
- `--L`: Loss Cost.
- `--I`: ILS Cost.
- `--V`: Verbose Mode[1=True].

## Output

The tool generates a CSV output file with the following columns:

- Process
- Replicate
- Labeled Gene_tree
- Species_Tree
- Number of Duplications
- Number of NNI (Nearest Neighbor Interchange)
- Number of Losses

Example entry in the output file:

```
Process,   Replicate,  Gene_tree,   Species_Tree,    Duplication,   NNI,   Loss
reconcILS,     0,     "(B,(C,A)1I  C);", "(A,(B,C));"    ,0,             1,     0
```

## Example Usage

You can use reconcILS as follows:

```bash
python ./reconcILS/reconcILS.py --spTree '(A,(B,C));' --gTree '((A,C),B);' --output 'result.csv'
```

