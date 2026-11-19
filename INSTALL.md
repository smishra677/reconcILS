# Install reconcILS

## Requirements

- Python 3.9 or later
- pandas
- matplotlib
- python-igraph
- Biopython
- ETE3
- Rich

## Recommended: virtual environment

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install .
```

On Windows, activate the environment with `.venv\Scripts\activate`.

## Verify the command

```bash
reconcils --help
```

Run the included example:

```bash
reconcils \
  --sp_file examples/spTree.tre \
  --gt examples/gTree.tre \
  --output results.csv \
  --verbose 1
```

Species-tree leaves must be `SpeciesName`; gene-tree leaves must be `SpeciesName_copy`. The run creates a deterministic two-character internal map and restores the original labels in its outputs.

## Editable installation

For development from a source checkout:

```bash
python -m pip install -e .
```

Inspect the parameterized parallel runner with:

```bash
python examples/parallel_reconciliation.py --help
```

## Uninstall

```bash
python -m pip uninstall reconcils
```

## Conda installation

See [CONDA.md](CONDA.md) for environment creation and conda-build packaging.
