# reconcILS: Simulation and Reconciliation Analyses ( rearranged by AI)
Please find Experiment_V6 in dryad: 


This directory contains the commands and scripts used for the simulation and reconciliation analyses in **reconcILS: A gene tree-species tree reconciliation algorithm that allows for incomplete lineage sorting**.

reconcILS is a parsimony algorithm for gene tree-species tree reconciliation that accounts for **duplication, loss, and incomplete lineage sorting (ILS)**. It uses **nearest-neighbor interchange (NNI)** rearrangements to model discordance due to the coalescent process. The analyses compare reconcILS with duplication-loss reconciliation in **ETE3** and with **DLCpar**. Gene trees were simulated using **dupcoal** under the **multispecies coalescent with duplication and loss (MSC-DL)** model and using **SimPhy**.

## Directory structure

The experimental scripts and input files are expected inside the `experiments/` directory of the `reconcILS` repository. The commands in `generate.sh` use relative script and tree paths and should be run from that directory:

```text
reconcILS/
└── experiments/
    ├── README.md
    ├── generate.sh
    ├── generate.r
    ├── scale.py
    ├── simulate_noise.py
    ├── test_reconcILS_default_40.py
    ├── trigger_dlcpar.py
    ├── trigger_simphy.py
    ├── run_all_paralogy.py
    ├── species_40_with_branch_lengths.tre
    ├── 3sp.tre
    └── Data_with_duplication_40_species_0.03_0.03/
        ├── trees.tre
        └── summary.csv
```

The tree files and dataset directory shown above are inputs to the corresponding experiments; they must be provided or generated before running those commands. Other datasets and species-tree files are required for the additional analyses described below.

```bash
cd reconcILS/experiments
```

The original `generate.sh` includes absolute paths under `/N/slate/.../reconcILS/experiments/`, including paths to dupcoal and DLCpar installations. These paths refer to the computing environment used for the analyses and must be replaced with the corresponding paths on another system. The `generate.r` output paths must also be adjusted before running it elsewhere.

## Scripts

| Script | Analysis |
| --- | --- |
| `generate.sh` | Commands for simulation and reconciliation analyses |
| `generate.r` | Generate 10 random species trees with 40 tips |
| `scale.py` | Scale species-tree branch lengths |
| `test_reconcILS_default_40.py` | Reconcile simulated gene trees using reconcILS and ETE3 |
| `simulate_noise.py` | Add gene tree error using NNI rearrangements |
| `trigger_dlcpar.py` | Compare reconcILS, DLCpar, and ETE3 |
| `trigger_simphy.py` | Reconcile gene trees simulated using SimPhy |
| `run_all_paralogy.py` | Analyze the primate gene trees |

## Simulation using dupcoal

The main analyses use **1,000 gene trees simulated on a 40-tip species tree** with duplication, loss, and coalescence. The duplication and loss rate parameters in `generate.sh` are both `0.03`.

```bash
python3 /path/to/dupcoal_v3.0.0.py \
  --stree species_40_with_branch_lengths.tre \
  --mu_par 0.03 --lambda_par 0.03 \
  --reps 1000 \
  --output Data_with_duplication_40_species_0.03_0.03
```

Reconcile the simulated gene trees using the default event costs: **duplication = 1.1, loss = 1, NNI = 1**.

```bash
python3 test_reconcILS_default_40.py \
  --data Data_with_duplication_40_species_0.03_0.03 \
  --species-tree species_40_with_branch_lengths.tre \
  --duplication-cost 1.1 --ils-cost 1 --loss-cost 1 \
  --workers 20
```

The input dataset contains `trees.tre` and `summary.csv`. The analyses compare inferred and simulated numbers of duplications, losses, and NNIs (**Figure 3; Supplementary Table S1**).

## Incomplete lineage sorting without duplication or loss

To examine discordance due solely to ILS, the manuscript uses simulations with no duplication or loss on species trees with 3 and 40 tips. The three-species command is:

```bash
python3 /path/to/dupcoal_v3.0.0.py \
  --stree 3sp.tre --mu_par 0 --lambda_par 0 \
  --reps 1000 \
  --output Data_with_duplication_3_species_0_0
```

The 40-species ILS-only dataset is named `Data_with_duplication_40_species_just_ILS` in `generate.sh`.

## Duplication and loss rates and event costs

The commands include simulations with unequal duplication and loss rate parameters:

| `mu_par` | `lambda_par` |
| ---: | ---: |
| 0.03 | 0.01 |
| 0.01 | 0.03 |
| 0.0001 | 0.03 |
| 0.03 | 0.0001 |

The reconciliation analyses also examine the effects of different event costs:

| Duplication | NNI | Loss |
| ---: | ---: | ---: |
| 1.1 | 1 | 1 |
| 3 | 1 | 4 |
| 4 | 1 | 3 |
| 1 | 5 | 1 |

These analyses are described in the **Appendix**.

## Long species-tree branches

To examine conditions with very little ILS, the manuscript uses species-tree branch lengths multiplied by **1,000**, with duplication and loss rate parameters of `0.00003`.

```bash
python3 scale.py \
  --input species_40_with_branch_lengths.tre \
  --scale 1000 \
  --output species_40_with_branch_lengths_1000.tre

python3 /path/to/dupcoal_v3.0.0.py \
  --stree species_40_with_branch_lengths_1000.tre \
  --mu_par 0.00003 --lambda_par 0.00003 \
  --reps 1000 \
  --output Data_with_duplication_40_species_0.00003_0.00003_x1000
```

This analysis is reported in **Supplementary Table S1**. An earlier command in `generate.sh` uses `--scale 100`; the manuscript describes the long-branch analysis as **x1000**.

## Ten additional species trees

The manuscript repeats analyses on **10 additional, randomly drawn species trees of varying heights** (**Supplementary Table S2**).

`generate.r` uses `ape::rtree()` to generate 10 rooted trees with 40 tips, assigns branch lengths using Grafen's method, and multiplies the branch lengths of each tree by a randomly sampled factor from 3 to 50. The trees are written to `species_40_scaled_10_simphy.tre` and the factors to `scale_factors_simphy.txt`.

```bash
Rscript generate.r
```

`generate.sh` contains separate loops for the 10 dupcoal datasets and the 10 SimPhy datasets. The former uses `species_40_scaled_10.tre`; the latter uses `species_40_scaled_10_simphy.tre`.

## Gene tree error

To assess the effects of gene tree estimation error, the analyses add **5, 10, or 15 errors per tree, on average**, using NNI rearrangements (**Supplementary Table S1**).

```bash
python3 simulate_noise.py --file Data_with_duplication_40_species_0.03_0.03/ --error 5
python3 simulate_noise.py --file Data_with_duplication_40_species_0.03_0.03/ --error 10
python3 simulate_noise.py --file Data_with_duplication_40_species_0.03_0.03/ --error 15
```

The resulting gene trees are named `error_trees_5.tre`, `error_trees_10.tre`, and `error_trees_15.tre`. They are reconciled using `test_reconcILS_default_40.py` with the corresponding `--file_name` argument.

## Comparison with DLCpar

The **Appendix** compares reconcILS with DLCpar using the three-species and fungi simulations. `generate.sh` includes the following three-species analysis:

```bash
python3 trigger_dlcpar.py \
  --species-tree 3sp.tre \
  --duplication-cost 1.1 --ils-cost 1 --loss-cost 1 \
  --workers 20 \
  --data Data_with_duplication_3_0.15_0.3 \
  --dlcpar-python /path/to/python2
```

The script also includes the three-species dataset with rate parameters `0.1` and `0.3`, and the fungi dataset with rate parameters `0.0015` and `0.003`.

## Comparison with SimPhy

Gene trees simulated using **SimPhy** are analyzed using `trigger_simphy.py`. The script compares duplication and loss inferences from reconcILS and ETE3. SimPhy does not report the number of NNIs per gene tree.

```bash
python3 trigger_simphy.py \
  --data /path/to/simphy_dataset \
  --workers 20 \
  --out /path/to/simphy_results
```

The 10-species-tree SimPhy analyses use the trees in `species_40_scaled_10_simphy.tre`.

## Analysis of primate genomes

The manuscript analyzes gene trees from **23 primate genomes**: **1,820 single-copy orthologs** and **11,555 gene trees containing both orthologs and paralogs** (**Figure 4**).

```bash
python3 run_all_paralogy.py \
  --genes ZeroCol_ASTRAL_ML_SCO_MIN27.tre \
  --mode sco \
  --species sp_tree_pruned.tre \
  --output-dir sco_efficient
```

The full set of gene trees can be analyzed using `--mode paralogy` and the corresponding input gene-tree file.

## Event labels and results

reconcILS reports **duplication**, **loss**, and **NNI** counts. The event labels are `D` (duplication), `L` (loss), and **`N` (NNI)**. The output column for NNI events is named `NNI`; the `--ils-cost` option specifies the NNI cost.

The experiment scripts produce reconciliation results, labeled species trees, and timing and memory measurements. The simulation analyses compare inferred events with the events recorded by dupcoal or SimPhy.

## Citation

Mishra, S. R., Smith, M. L., and Hahn, M. W. **reconcILS: A gene tree-species tree reconciliation algorithm that allows for incomplete lineage sorting.** bioRxiv. https://doi.org/10.1101/2023.11.03.565544

**dupcoal:** https://github.com/meganlsmith/dupcoal
