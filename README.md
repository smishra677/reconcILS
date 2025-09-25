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

Install the package via Conda:

```bash
conda install -c smishra677 reconcils
```




## Arguments 
| Argument       | Description                                       | Required | Default Value |
| -------------- | ---------------------------------                 | -------- | ------------- |
| `--spTree`     | Species tree in string // Location to Species Tree                        | Yes      | N/A           |
| `--gTree`      | Gene tree in string // Location to Gene Tree                           | Yes      | N/A           |
| `--output`     | Name of the output file with the extension .csv | Yes      | N/A           |
| `--D`          | Duplication Cost                                 | No       | 1.1           |
| `--L`          | Loss Cost                                       | No       | 1.0           |
| `--I`          | ILS Cost                                        | No       | 1.0           |
| `--V`          | Verbose Mode                                    | No       | 0             |
| `--F`          | Input as file(--F 1 for file input)                                   | No       | 0            |

## Output

The tool generates a CSV output file with the following columns:

- Process
- Replicate (By Default :0)
- Labeled Gene_tree
- Species_Tree
- Number of Duplications
- Number of NNI (Nearest Neighbor Interchange)
- Number of Losses

Example entry in the output file:

output.csv
| Process    | Replicate | Gene_tree      | Species_Tree   | Duplication | NNI | Loss |
|------------|-----------|----------------|----------------|-------------|-----|------|
| reconcILS  | 0         | "((C,A)1I C,B);	" | "((0-0-0   C,0-0-0   B)0-1-0   ,0-0-0   A)0-0-0   ;"   | 0           | 1   | 0    |

output_log.csv

| Gene_Tree     | Species_Tree  | Duplication_cost | NNI_cost | Loss_cost |
|-------------- | ------------- | ---------------- | --------  | --------- |
| "(B,(C,A));"  | "(A,(B,C));" | 1.1              | 1.0     | 1         |



## Parsing Output

### Labeled Species Tree
((0-0-0 C, 0-0-0 B) 0-1-0, 0-0-0 A) 0-0-0

Here, reconcILS produces a labeled species tree where each branch is labeled with the number of **Duplication-NNI-Loss** events for that branch. For instance, `((0-0-0 C, 0-0-0 B) 0-1-0, 0-0-0 A) 0-0-0` indicates there was one NNI move on the branch leading from `(B, C)` to `(A, (B, C))` on the species tree.

### Labeled Gene Tree
((C, A) 1I C, B)

Here, reconcILS labels the gene tree with Duplication (D), NNI (I), and Loss (L) events. For example, `((C, A) 1I C, B)` indicates there was one NNI move on `C` at the branch `(C, A)` to `(A, (B, C))` on the gene tree.

These labeled trees can be visualized with any Newick visualizer.


## Experiments:

For experiments, please see the experiment branch. (https://github.com/smishra677/reconcILS/blob/Experiments/reconcILS/Experiments.md)

## Example Usage

You can use reconcILS as follows:

Input as String:
```bash
python ./reconcILS/reconcILS.py --spTree '(A,(B,C));' --gTree '((A,C),B);' --output 'result.csv'
```

Input as file:
```bash
python ./reconcILS/reconcILS.py --spTree './spTree.tre' --gTree './gTree.tre' --output 'result.csv' --F 1
```
Using reconcILS from a Python file:

Please see the example.md at https://github.com/smishra677/reconcILS/blob/main/reconcILS/example/example.md


## Example Usage Conda


Use the package programmatically from Python:

```python

import sys
from reconcILS import reconcILS

# Simulate command-line arguments
sys.argv = [
    "reconcILS.py",
    "--spTree", "(A,(B,C));",     # Species tree in Newick format
    "--gTree", "((A,C),B);",      # Gene tree in Newick format
    "--output", "result.csv"      # Output file name
]

# Run the reconciliation
reconcILS.main()
```

This will generate:
- `result.csv`: Reconciliation summary including duplications, losses, and ILS events
- `result_log.csv`: Metadata log including input trees and cost parameters
- `resultspecies_acr.json` and `resultacr_species.json`: Mappings of species names to acronyms (if needed)

> This interface allows calling `reconcILS` directly from Python scripts or notebooks (e.g., Google Colab), without using the command line.
