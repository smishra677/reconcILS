# reconcILS — User Guide

**reconcILS** is a Python tool for reconciling gene trees with a species tree and reporting duplication, incomplete lineage sorting (ILS/NNI), and loss events. It supports a single gene tree through the `reconcils` command and multiple gene trees through a parameterized parallel runner.

This guide is intended for **new users**: it explains how to install the software, prepare your own data, run an analysis, interpret the generated files, and troubleshoot common problems. It follows the practical user-guide organization of DAFT.

## What reconcILS does

```text
Rooted species tree (SpeciesName) + rooted gene tree(s) (SpeciesName_copy)
                             |
                             v
                   Parse and validate inputs
                             |
                             v
              Deterministic internal species aliases
                             |
                             v
                 reconcILS reconciliation engine
                             |
                             v
            Event counts, labeled trees, logs, alias map
                             |
                             v
                Restore original species/copy labels
```

The reconciliation engine uses event costs for duplication (`D`), loss (`L`), and ILS (`I`). The result CSV uses the column `NNI` for ILS-related events. Differences between gene-tree and species-tree topology are **not** treated as input-validation errors.

## 1. Get the source code

Download the repository ZIP from its GitHub page and extract it, or clone the repository using its GitHub **Code → HTTPS** URL:

```bash
git clone <repository-HTTPS-URL>
cd <cloned-repository-directory>
```

The project includes `reconcILS/` (Python package), `examples/` (sample trees and parallel runner), `tests/`, `environment.yml`, `conda-recipe/`, and installation documentation. The GitHub repository URL is intentionally not hardcoded here; use the URL of the repository where this source is published.

## 2. Install reconcILS

### Recommended: Conda

Install Conda or Miniforge first, then run from the repository root:

```bash
conda env create -f environment.yml
conda activate reconcils
reconcils --help
```

### Alternative: pip / virtual environment

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
reconcils --help
```

Dependencies include pandas, matplotlib, igraph, Biopython, ETE3, and Rich. For more detail see [INSTALL.md](INSTALL.md) and [CONDA.md](CONDA.md).

## 3. Prepare input trees

**Species-tree leaves** must contain species names, with no copy suffix. **Gene-tree leaves** must be `SpeciesName_copy`, with a nonempty copy identifier. Every full gene-copy label must be unique *within a gene tree*.

Species tree (`species.tre`):

```newick
((Homo_sapiens,Pan_troglodytes),Mus_musculus);
```

Gene tree (`gene.tre`):

```newick
((Homo_sapiens_copy1,Mus_musculus_copy1),Pan_troglodytes_copy1);
```

Multiple gene trees (`gene_trees.tre`), one per line:

```newick
((Homo_sapiens_copy1,Mus_musculus_copy1),Pan_troglodytes_copy1);
((Homo_sapiens_copy2,Pan_troglodytes_copy2),Mus_musculus_copy2);
```

Species names can contain underscores. reconcILS matches the **longest exact species-tree name followed by `_`** and retains the remainder as the copy identifier. For example, `Homo_sapiens_copy_17` is interpreted as species `Homo_sapiens` and copy `copy_17`.

Use complete, semicolon-terminated Newick trees. The validator expects binary trees and unique species names. Branch lengths and internal node names are not required for the topology-based workflow. Gene trees may contain fewer species than the species tree, but may not contain unknown species. A gene tree can contain multiple distinct copies of the same species.

### Automatic internal aliases

Species names are sorted lexicographically, then assigned deterministic two-character identifiers such as `A0`, `A1`, `A2`. A gene leaf `Homo_sapiens_copy1` becomes `A0_copy1` *inside the algorithm* and is restored in user-facing results. Parsed tree leaves are renamed structurally; the program does not use substring replacement to translate names.

A JSON mapping report is saved for reproducibility. Advanced users can provide `--mapping` with JSON, CSV, or TSV (`species,alias` columns). Internal aliases must be unique, two-character uppercase letters/digits, without underscores. The namespace supports at most 1,296 distinct species.

## 4. Quick start: one gene tree

From the repository root:

```bash
reconcils \
  --sp_file examples/spTree.tre \
  --gt examples/gTree.tre \
  --output results.csv \
  --verbose 1
```

You can also provide a species tree inline:

```bash
reconcils \
  --sp '((Homo_sapiens,Pan_troglodytes),Mus_musculus);' \
  --gt '((Homo_sapiens_copy1,Mus_musculus_copy1),Pan_troglodytes_copy1);' \
  --output results.csv
```

The CLI detects whether `--gt` is an existing file or inline Newick. **Do not pass the legacy `--F` option.** Run `reconcils --help` to see the installed command options.

### Single-tree arguments

| Option | Meaning | Default |
|---|---|---|
| `--sp` / `--sp_file` | Species tree inline or in a file (choose one) | Required |
| `--gt` | Gene tree inline or file | Required |
| `--output` | Result CSV filename | Required |
| `--D` | Duplication cost | 1.1 |
| `--L` | Loss cost | 1 |
| `--I` | ILS cost | 1 |
| `--ils-timeout` | Internal ILS-search timeout, seconds | 5000 |
| `--verbose` | `0` plain output; `1` Rich display | 0 |
| `--mapping` | Optional alias map | Generated automatically |
| `--alias-output` | Override alias JSON output | Derived from result filename |
| `--species-format`, `--gene-format` | ETE3 input format integers | 1 |

## 5. Run many gene trees in parallel

The runner in `examples/parallel_reconciliation.py` takes **your own** species tree and gene-tree file. No simulation data, predefined taxa, or `summary.csv` are required.

```bash
python examples/parallel_reconciliation.py \
  --species-tree species.tre \
  --gene-trees gene_trees.tre \
  --output batch_results.csv \
  --workers 4 \
  --verbose 1
```

A more explicit run:

```bash
python examples/parallel_reconciliation.py \
  --species-tree species.tre \
  --gene-trees gene_trees.tre \
  --output batch_results.csv \
  --duplication-cost 1.1 \
  --loss-cost 1 \
  --ils-cost 1 \
  --workers 8 \
  --library-threads 1 \
  --ils-timeout 600 \
  --timeout 1800 \
  --checkpoint-every 100 \
  --verbose 1
```

`--workers` controls independent worker processes. `--library-threads` limits native numerical-library threads per worker. `--limit` restricts the number of input gene trees. The parallel runner requires Linux or macOS for its POSIX file-locking mechanism.

## 6. Timeouts: two separate controls

| Option | Where available | What it limits | Default |
|---|---|---|---|
| `--ils-timeout SECONDS` | Single-tree CLI and parallel runner | One internal ILS search inside the reconciliation algorithm | 5000 seconds |
| `--timeout SECONDS` | Parallel runner | Overall wall-clock processing time for one gene tree, including worker startup | No outer limit |

Example: `--ils-timeout 600 --timeout 1800` requests a 10-minute internal ILS search limit and a 30-minute overall per-tree limit. The **outer** timeout kills a worker that exceeds its limit; the **internal** timeout uses the algorithm's existing `ThreadPoolExecutor` timeout, which does not forcibly terminate a running thread. For an enforceable overall bound on a difficult tree, use the parallel runner's `--timeout`.

Both timeout values must be positive finite numbers. Timeout settings are included in the parallel checkpoint metadata, so resuming with a different internal ILS timeout requires a new state/output directory. `--timeout` is a scheduling limit and may be changed when rerunning.

## 7. Validation and logs

Before reconciliation, the tool checks Newick syntax, semicolon termination, binary topology, duplicate species names, duplicate full gene-copy labels, copy suffixes, and whether each gene leaf maps to a species in the supplied species tree. The validator **does not reject** gene/species topology discordance.

A failed single-tree run reports the input source and reason. The parallel runner identifies the gene-tree index and retains failure diagnostics under its state directory. Use `--verbose 1` for Rich panels and more algorithm diagnostics; `--worker-logs` keeps logs for successful parallel trees too.

## 8. Output files and how to read them

For `reconcils --output results.csv`, expect:

| File | Purpose |
|---|---|
| `results.csv` | Reconciliation row with original gene/species names and event counts |
| `results_species_aliases.json` | Forward and reverse deterministic species/alias mapping |
| `results.csv_validation.json` | Input-validation summary |
| `results.csv_taxon_report.tsv` | Taxon and gene-copy report |
| `results.csv_reconcILS_log.txt` | Validation/run log |
| `results_log.csv` | Input and event-cost log |

The main result CSV includes `Gene_tree`, `Species_Tree`, `Duplication`, `NNI`, and `Loss`. `Duplication` and `Loss` count their respective reconciliation events; `NNI` is the reported ILS/NNI event count. A count is an algorithm output, not an input-validation judgment. Original species and copy names are restored before writing user-facing results.

For the parallel runner, `--output` controls the result CSV; it also writes an aggregate labeled species tree and a species alias JSON. Use `--labeled-species-output`, `--alias-output`, and `--state-dir` to control their locations. Parallel CSV rows include a `Replicate` index and `Seconds` for runtime.

## 9. Checkpointing and resuming parallel runs

Parallel processing writes a journal and checkpoint files to `--state-dir` (by default, a hidden directory beside the output). Rerun the **same command** to resume compatible work. The journal checks the input trees, species mapping, costs, internal ILS timeout, formats, and algorithm fingerprint. If any of these change, choose a new output/state directory. Use `--retry-failed` to rerun failed records. The result CSV is written in input order; a failure may prevent later completed rows from appearing until that failure is resolved.

## 10. Rich terminal display

Add `--verbose 1` to the single-tree CLI or parallel runner. The single-tree command shows input/validation panels and event summaries. The parallel runner shows progress, active workers, saved trees, failures, and checkpoint status. Rich is a display layer; it does not change the reconciliation algorithm.

## 11. Troubleshooting

| Symptom | Suggested action |
|---|---|
| `reconcils: command not found` | Activate the environment; verify with `python -m pip show reconcils` and `which reconcils` |
| `No module named ete3` or `igraph` | Install the dependencies in the same environment used to run reconcILS |
| `Gene leaf ... is missing a copy ID` | Rename it from `SpeciesName` to `SpeciesName_copy1` (or another unique suffix) |
| Unknown gene species | Check that the gene leaf begins with an exact species-tree leaf name followed by `_` |
| Duplicate copy label | Give each gene copy a unique full label within its gene tree |
| Empty event records / internal timeout | Inspect worker logs, check the inputs and costs, and consider `--ils-timeout` / `--timeout` |
| Journal metadata mismatch | Use a fresh `--state-dir` and output path when changing inputs, costs, or internal timeout |
| No Rich dashboard | Run with `--verbose 1` in an interactive terminal |

## 12. Suggested reporting and reproducibility

When reporting reconcILS results, retain the species and gene tree inputs, version of reconcILS, duplication/loss/ILS costs, timeout settings, alias JSON, validation report, and result CSV. For parallel analyses, retain the checkpoint metadata and record the number of successfully processed trees. This makes it possible to reproduce how gene copies were matched to species and which reconciliation settings were used.

For detailed examples see [examples/README.md](examples/README.md); for installation and Conda packaging see [INSTALL.md](INSTALL.md) and [CONDA.md](CONDA.md).

## Output organization (DAFT style)

For `--output analysis.csv`, reconcILS automatically organizes generated files
into three sibling directories beside the requested output path:

```text
reconcILS_results_analysis/   # result CSV and labeled species tree (parallel)
reconcILS_extras_analysis/    # species alias map and validation reports
reconcILS_log_analysis/       # run logs and parallel checkpoint state
```

The single-tree command organizes files after a successful run. The parallel
runner writes directly into these folders so its checkpoint journal remains
stable across interruptions and retries. Re-run with the same output name to
resume parallel processing; use a different output name for a separate run.
No unrelated files are moved.
