# Conda installation and packaging

The project includes a development environment and a `conda-build` recipe. Runtime dependencies are resolved from `conda-forge`, including `python-igraph`, Biopython, ETE3, pandas, matplotlib, and Rich.

## Create a development environment

```bash
conda env create -f environment.yml
conda activate reconcils
reconcils --help
```

Run inline trees using the required species/copy convention:

```bash
reconcils \
  --sp '((Homo_sapiens,Pan_troglodytes),Mus_musculus);' \
  --gt '((Homo_sapiens_copy1,Mus_musculus_copy1),Pan_troglodytes_copy1);' \
  --output result.csv \
  --verbose 1
```

Or use the included files:

```bash
reconcils \
  --sp_file examples/spTree.tre \
  --gt examples/gTree.tre \
  --output result.csv
```

The wrapper detects files automatically; users do not pass the legacy `--F` option. It assigns deterministic two-character internal IDs and writes original `SpeciesName` and `SpeciesName_copy` labels in output files.

## Run the parameterized parallel example

```bash
python examples/parallel_reconciliation.py \
  --species-tree /path/to/species.tre \
  --gene-trees /path/to/gene_trees.tre \
  --output /path/to/reconcILS_result.csv \
  --workers 4 \
  --verbose 1
```

See [examples/README.md](examples/README.md) for mapping, checkpoint, timeout, cost, and output parameters.

## Build and install a Conda package

```bash
conda env create -f environment-build.yml
conda activate reconcils-build
conda build conda-recipe -c conda-forge
```

Install the local package in a clean environment:

```bash
conda create -n reconcils-test -c conda-forge python=3.11
conda activate reconcils-test
conda install -c conda-forge --use-local reconcils
reconcils --help
```

If Conda cannot find the local package, install the explicit package path printed by `conda build`.

## Notes

- `python-igraph` supplies the `igraph` Python module.
- `fcntl` is available on POSIX systems rather than native Windows; the parallel runner uses it for output locking.
- The source distribution also supports `python -m pip install .`.
- The recipe is `noarch: python`; Conda resolves platform-specific dependencies when installed.
- Verify `conda build` on the target platform before publishing the package to a channel.
