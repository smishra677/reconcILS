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
| `--spTree`     | Species tree in string                           | Yes      | N/A           |
| `--gTree`      | Gene tree in string                              | Yes      | N/A           |
| `--output`     | Name of the output file with the extension .csv | Yes      | N/A           |
| `--D`          | Duplication Cost                                 | No       | 1.1           |
| `--L`          | Loss Cost                                       | No       | 1.0           |
| `--I`          | ILS Cost                                        | No       | 1.0           |
| `--V`          | Verbose Mode                                    | No       | 0             |

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
| Process    | Replicate | Gene_tree      | Species_Tree   | Duplication | NNI | Loss |
|------------|-----------|----------------|----------------|-------------|-----|------|
| reconcILS  | 0         | "(B,(C,A)1I  C);" | "(A,(B,C));"   | 0           | 1   | 0    |


## Example Usage

You can use reconcILS as follows:

```bash
python ./reconcILS/reconcILS.py --spTree '(A,(B,C));' --gTree '((A,C),B);' --output 'result.csv'
```

