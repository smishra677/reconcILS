# Experiments V6

Scripts, species trees, and results for the [reconcILS paper](https://doi.org/10.1101/2023.11.03.565544). Source code: [reconcILS](https://github.com/smishra677/reconcILS); simulator: [dupcoal](https://github.com/meganlsmith/dupcoal).

## Dataset

### Input

1. Clone `reconcILS` and extract this folder directly into the repository as `reconcILS/Experiments_V6/`, alongside `reconcILS.py` and `utils_reconcILS/`.
2. [all_species tree/](all_species%20tree/) contains the three- and 40-species trees, the ten trees used for Table S2, and the pruned primate species tree.
3. [smith_et_al_2022_primates/](smith_et_al_2022_primates/) contains the single-copy ortholog (SCO) and all-paralogy gene trees and their reconciliation results.

Run the commands below from `reconcILS/Experiments_V6/`. The Python scripts import reconcILS from the parent directory; keep generated data folders inside `Experiments_V6/`.

The original simulated gene-tree datasets are not included but can be provided upon request.. Generate them with dupcoal or SimPhy before running the simulation analyses.

## Experiments

### Requirements

- Install `reconcILS` and its dependencies, including `ete3`.
- Install [dupcoal](https://github.com/meganlsmith/dupcoal) to simulate gene trees.
- Install `DLCpar` for DLCpar comparisons.
- Install `SimPhy` for SimPhy comparisons.
- Install R and the `ape` package to run `generate.r`.

### Simulated dataset

A 40-species example is:

```bash
python3 /path/to/dupcoal_v3.0.0.py \
  --stree "all_species tree/species_40_with_branch_lengths.tre" \
  --mu_par 0.03 --lambda_par 0.03 --reps 1000 \
  --output Data_with_duplication_40_species_0.03_0.03

python3 test_reconcILS_default_40.py \
  --data Data_with_duplication_40_species_0.03_0.03 \
  --species-tree "all_species tree/species_40_with_branch_lengths.tre"
```

The generated data folder must contain `trees.tre` and `summary.csv`. Other scripts are [trigger_dlcpar.py](trigger_dlcpar.py), [trigger_simphy.py](trigger_simphy.py), and [simulate_noise.py](simulate_noise.py) for DLCpar, SimPhy, and gene-tree error analyses.

### All-paralogy / single-copy ortholog (SCO) dataset

All-paralogy is the default. Run the following from `reconcILS/Experiments_V6/`:

```bash
python3 run_all_paralogy.py
```

For SCO, specify the SCO gene trees:

```bash
python3 run_all_paralogy.py \
  --genes smith_et_al_2022_primates/ZeroCol_ASTRAL_ML_SCO_MIN27.tre
```

The script selects SCO mode from the input filename. Results are written to `paralogy_efficient/` or `sco_efficient/` by default.

## Generation and analysis scripts

All experiment-generation and analysis commands are listed in [generate.sh](generate.sh), including simulations, cost settings, method comparisons, gene-tree error tests, and primate analyses.

- [generate.r](generate.r): Generates ten random 40-species trees and scaled versions for SimPhy.
- [scale.py](scale.py): Scales species-tree branch lengths for the long-branch experiments.
- [genrate_simphy.sh](genrate_simphy.sh): Generates SimPhy datasets.
- [test_reconcILS_default_40.py](test_reconcILS_default_40.py): Runs reconcILS and ETE3 on simulated gene trees.
- [trigger_dlcpar.py](trigger_dlcpar.py): Runs DLCpar comparisons.
- [trigger_simphy.py](trigger_simphy.py): Runs reconciliations on SimPhy datasets.
- [simulate_noise.py](simulate_noise.py): Adds NNI errors to simulated gene trees.
- [run_all_paralogy.py](run_all_paralogy.py): Runs primate analyses, using all-paralogy by default or SCO when specified.

## Outputs

- Simulations: `reconcILS_result_*.csv` contains the true event counts and reconcILS/ETE3 estimates. `reconcILS_time_result_*.csv`, `reconcILS_memory_result_with_dups_loss_*.csv`, and `reconcILS_labeled_sp_*.csv` contain timing, memory, and species-tree event labels.
- Primates: `bio_result_*_pruned67.csv` contains per-gene event counts; `labeled_L_*_pruned67.csv` contains species-tree event labels; `*_pruned_Gene_trees67.csv` contains the processed gene trees. Timing, memory, and status files are also saved. Existing outputs are in [smith_et_al_2022_primates/](smith_et_al_2022_primates/).

## Figures and other files

| Paper | Files |
| --- | --- |
| Figure 3; Figure S11 | [40-species script](test_reconcILS_default_40.py), [notebook](Produce_results-40_species_V_6.ipynb) |
| Figure S10 | [3-species tree](all_species%20tree/3sp.tre), [historical simulation commands](generate.sh) |
| Figure 4; Figure S12 | [Primate results](smith_et_al_2022_primates/), [analysis script](run_all_paralogy.py) |
| Figure S9; Table S1 | [40-species trees](all_species%20tree/), [error script](simulate_noise.py), [Table S1](paper/supplementary_table_S1.pdf) |
| Table S2 | [Ten species trees](all_species%20tree/species_40_scaled_10.tre), [SimPhy-scaled trees](all_species%20tree/species_40_scaled_10_simphy.tre), [Table S2](paper/supplementary_table_S2.pdf) |

The final ten-tree SimPhy commands in [generate.sh](generate.sh) use `species_40_scaled_10_simphy.tre`. The bundled paper trees are in [all_species tree/](all_species%20tree/); running `generate.r` creates new trees, not the exact paper inputs. The notebooks also contain original computing-environment paths. See the [paper PDF](paper/reconcILS_preprint_v6.pdf) for details.
