# reconcILS (for old reconcILS main/reconcILS-old)
for experiments related to paper:  https://github.com/smishra677/reconcILS/tree/Experiments

**Gene-tree reconciliation with duplication, incomplete lineage sorting, and gene loss**

reconcILS reconciles rooted gene trees against a rooted species tree. It reports inferred duplication, incomplete lineage sorting (ILS/NNI), and loss events, and provides a command-line interface for individual trees and a parallel workflow for larger collections.

## Contents

- [Installation](#installation)
- [Input format](#input-format)
- [Quick start](#quick-start)
- [Parallel analysis](#parallel-analysis)
- [Parameters and timeouts](#parameters-and-timeouts)
- [Output files](#output-files)
- [Species and gene-copy labels](#species-and-gene-copy-labels)
- [Validation and troubleshooting](#validation-and-troubleshooting)
- [Reproducibility](#reproducibility)
- [License](#license)

## Installation

Python 3.9 or newer is required. Clone the repository and install from its root directory.

**Conda (recommended):**

```bash
git clone https://github.com/smishra677/reconcILS.git
cd reconcILS
conda env create -f environment.yml
conda activate reconcils
reconcils --help
```

**Python virtual environment:**

```bash
git clone https://github.com/smishra677/reconcILS.git
cd reconcILS
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
reconcils --help
```

The package uses ETE3, Biopython, python-igraph, pandas, matplotlib, and Rich. For installation details, see [INSTALL.md](INSTALL.md) and [CONDA.md](CONDA.md). The Conda environment installs this checkout; it does not automatically follow later repository commits.

## Input format

Inputs are rooted, semicolon-terminated Newick trees. Species-tree leaves contain **species names**, while gene-tree leaves contain **species names followed by an underscore and a copy identifier**.

**Species tree (`species.tre`):**

```newick
((Homo_sapiens,Pan_troglodytes),Mus_musculus);
```

**Gene tree (`gene.tre`):**

```newick
((Homo_sapiens_copy1,Mus_musculus_copy1),Pan_troglodytes_copy1);
```

A file of multiple gene trees (`genes.tre`) contains one Newick tree per line:

```newick
((Homo_sapiens_copy1,Mus_musculus_copy1),Pan_troglodytes_copy1);
((Homo_sapiens_copy2,Pan_troglodytes_copy2),Mus_musculus_copy2);
```

Species names may contain underscores. Gene leaves are matched against the longest complete species-tree name followed by `_`; the remaining suffix is retained as the copy identifier. Thus `Homo_sapiens_copy_17` belongs to `Homo_sapiens` and has copy identifier `copy_17`.

Species names must be unique within the species tree, and complete gene-copy labels must be unique within each gene tree. Gene trees may omit species, but cannot contain species absent from the species tree. The input validator checks binary topology; branch lengths and internal node labels are not needed for topology-based reconciliation.

## Quick start

Reconcile one gene tree using the example input files:

```bash
reconcils \
  --sp_file examples/spTree.tre \
  --gt examples/gTree.tre \
  --output analysis.csv \
  --verbose 1
```

Trees can also be provided directly:

```bash
reconcils \
  --sp '((Homo_sapiens,Pan_troglodytes),Mus_musculus);' \
  --gt '((Homo_sapiens_copy1,Mus_musculus_copy1),Pan_troglodytes_copy1);' \
  --output analysis.csv
```

The `--gt` option accepts either inline Newick or an existing file path. A separate file-mode flag is not needed. Use `reconcils --help` for the complete command reference.

## Parallel analysis

For multiple gene trees, use the parameterized runner. No simulated dataset or ground-truth summary file is required.

```bash
python examples/parallel_reconciliation.py \
  --species-tree species.tre \
  --gene-trees genes.tre \
  --output batch.csv \
  --workers 4 \
  --verbose 1
```

To configure event costs and execution limits:

```bash
python examples/parallel_reconciliation.py \
  --species-tree species.tre \
  --gene-trees genes.tre \
  --output batch.csv \
  --duplication-cost 1.1 \
  --loss-cost 1 \
  --ils-cost 1 \
  --workers 8 \
  --ils-timeout 600 \
  --timeout 1800 \
  --checkpoint-every 100
```

The runner uses independent worker processes, writes checkpoints, and can resume compatible runs. The parallel runner requires a POSIX environment (Linux or macOS) for file locking. See [examples/README.md](examples/README.md) for additional examples and options.

## Parameters and timeouts

### Single-tree command

| Option | Description | Default |
| --- | --- | --- |
| `--sp` / `--sp_file` | Species tree as Newick text or file (choose one) | Required |
| `--gt` | Gene tree as Newick text or file | Required |
| `--output` | Output CSV name | Required |
| `--D` | Duplication cost | 1.1 |
| `--L` | Loss cost | 1 |
| `--I` | ILS cost | 1 |
| `--ils-timeout` | Time limit for an internal ILS search, seconds | 5000 |
| `--mapping` | Optional species-to-alias mapping file | Generated |
| `--alias-output` | Mapping JSON filename | Derived from output |
| `--species-format`, `--gene-format` | ETE3 Newick format codes | 1 |
| `--verbose` | `0` for plain output, `1` for Rich terminal display | 0 |

### Additional parallel-runner options

| Option | Description | Default |
| --- | --- | --- |
| `--workers` | Number of worker processes | Up to 8, based on HPC(`SLURM_CPUS_PER_TASK`) |
| `--library-threads` | Native numerical-library threads per worker | 1 |
| `--limit` | Process only the first N trees | All |
| `--timeout` | Overall wall-clock time limit per tree, seconds | None |
| `--ils-timeout` | Internal ILS-search time limit, seconds | 5000 |
| `--checkpoint-every` | Number of completed records between exports | 100 |
| `--retry-failed` | Retry previously failed records | Off |
| `--worker-logs` | Retain logs for successful trees | Off |
| `--state-dir` | Checkpoint and journal directory | Generated |
| `--labeled-species-output` | Aggregated labeled species-tree output | Derived from output |

The two timeouts serve different purposes. `--ils-timeout` applies to an individual ILS search within the algorithm; its thread-based timeout does **not** forcibly stop the underlying thread. `--timeout` is available for parallel runs and terminates a worker process that exceeds the overall per-tree limit. Both values are in seconds and must be positive.

## Output files

Output files are organized by analysis name. For example, `--output analysis.csv` produces directories beside the requested output location:

```text
reconcILS_results_analysis/
    analysis.csv
reconcILS_extras_analysis/
    analysis_species_aliases.json
    analysis.csv_validation.json
    analysis.csv_taxon_report.tsv
reconcILS_log_analysis/
    analysis_log.csv
    analysis.csv_reconcILS_log.txt
```

The result CSV reports the input trees and inferred `Duplication`, `NNI`, and `Loss` counts. The `NNI` column reports the ILS/NNI events identified by reconciliation. The mapping JSON stores both directions of the species-name/alias mapping. The validation report and taxon report document the processed inputs; logs record run details.

Parallel runs additionally generate an aggregated labeled species tree and maintain journal/checkpoint data in the log directory. Paths can be customized with the parallel runner's output options. Existing files in destination folders are not silently overwritten by the single-tree organizer; use a new analysis name when appropriate.

## Species and gene-copy labels

Species names are assigned deterministic internal two-character identifiers, based on sorted original names. For example:

```text
Homo_sapiens       -> A0
Mus_musculus       -> A1
Pan_troglodytes    -> A2
```

A gene leaf such as `Homo_sapiens_copy1` becomes `A0_copy1` internally. Copy suffixes are preserved, and original names are restored in user-facing outputs. This mapping operates on parsed tree leaves rather than replacing substrings in Newick text.

The generated `*_species_aliases.json` file includes both `original_to_internal` and `internal_to_original` mappings. An optional mapping file may be supplied with `--mapping` (JSON, CSV, or TSV; tabular files use `species,alias` columns). Internal identifiers must be unique two-character uppercase alphanumeric values without underscores. The two-character namespace supports up to 1,296 species.

## Validation and troubleshooting

Validation runs automatically before reconciliation. It checks tree syntax, topology, duplicate labels, copy identifiers, and compatibility between gene-tree species and species-tree leaves. Gene/species topology discordance is not itself a validation error.

| Problem | Resolution |
| --- | --- |
| `reconcils: command not found` | Activate the installation environment and rerun `reconcils --help`. |
| Missing `ete3` or `igraph` | Install dependencies in the Python environment used to launch reconcILS. |
| Missing copy identifier | Use `SpeciesName_copy1`, not just `SpeciesName`, in the gene tree. |
| Unknown species | Match each gene leaf to an exact species-tree name followed by `_`. |
| Duplicate gene label | Assign a distinct full copy label to each leaf in a gene tree. |
| Internal ILS search timeout | Review logs and adjust `--ils-timeout` if appropriate. |
| Long-running parallel tree | Set `--timeout` for an enforceable worker-level limit. |
| Checkpoint metadata mismatch | Start a fresh output/state directory when changing inputs or algorithm settings. |

Use `--verbose 1` to display the Rich terminal interface. For parallel analyses, failed-tree diagnostics are retained with checkpoint state; `--worker-logs` also retains successful-tree logs.

## Reproducibility

For reproducible analyses, retain the original Newick files, reconcILS version, event costs, timeout settings, species-alias mapping JSON, and result CSV. Parallel analyses should also retain their checkpoint metadata and record the number of successfully processed gene trees. Re-running with different trees, mappings, costs, or internal ILS timeout settings may require a new checkpoint directory.

## License

reconcILS is distributed under the terms of the [MIT License](LICENSE).
