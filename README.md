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
| reconcILS  | 0         | "(B,(C,A)1I  C);" | "(A,(B,C));"   | 0           | 1   | 0    |

output_log.csv

| Gene_Tree     | Species_Tree  | Duplication_cost | NNI_cost | Loss_cost |
|-------------- | ------------- | ---------------- | --------  | --------- |
| "(B,(C,A));"  | "(A,(B,C));" | 1.1              | 1.0     | 1         |



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

