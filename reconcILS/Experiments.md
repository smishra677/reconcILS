# Dataset

## Input

### Folders: 10_30, zero_loss, zero_loss_zero_dups
- These folders contain gene trees simulated by DupCoal.
- The species tree used is inside the folder as "sp_tree.tre".
- The command given to DupCoal is also found inside the folder as "command.command".

### Files
- **ZeroCol_ASTRAL_ML_ALLPARALOGS_MIN27.tre, ZeroCol_ASTRAL_ML_SCO_MIN27.tre**: These are gene trees from Smith et al., 2020.
- **sp_tree.tre, sp.tre, sp_tree_pruned.tre**: These are the species trees used in our experiments.
  1. **sp_tree.tre**: Used in simulated experiments (10_30, zero_loss, zero_loss_zero_dups).
  2. **sp.tre**: The species tree from Smith et al., 2020.
  3. **sp_tree_pruned.tre**: Species tree from "sp.tre" after pruning outgroups and converting species into keys provided in the map.

# Experiments

## Requirements
1. Extract the experiment folder into `./reconcILS`.
2. Install `ete3` and `DLCpar`.

## Single Copy Orthologs / All Paralogy Dataset
3. Run the appropriate script: `python ./run_<SCO/all_paralogy>.py`.

## Simulated Dataset
4. Replace the `gene_folder` on line 107 with folder names from the input. For example: `./10_30`, `./zero_loss`, `./zero_loss_zero_dups`.

# Outputs

### Single Copy Orthologs / All Paralogy Dataset
By the end of the run, the program will output three types of files:

1. **Reconciliation Results**
   - `bio_result_ZeroCol_ASTRAL_ML_<SCO/all_paralogy>_MIN27_pruned.csv`
   - `labeled_L_ZeroCol_ASTRAL_ML_<SCO/all_paralogy>_MIN27_pruned.csv`

2. **Timing**
   - `timing_ZeroCol_ASTRAL_ML_<SCO/all_paralogy>_MIN27_pruned.csv`
   - `timing_ZeroCol_ASTRAL_ML_<SCO/all_paralogy>_MIN27_intro_1_write_pruned.csv`

3. **Output After Pruning Gene Trees**
   - `ZeroCol_ASTRAL_ML_<SCO/all_paralogy>_MIN27_pruned_Gene_trees.csv`

**Output Explanation:**
- `bio_result_ZeroCol_ASTRAL_ML_<SCO/all_paralogy>_MIN27_pruned.csv`: Records all the gene trees, with event numbers (i.e., number of duplications, NNI, and loss) required for reconciling each gene tree with the species tree for both ReconcILS and ete3.
- `labeled_L_ZeroCol_ASTRAL_ML_SCO_MIN27_pruned.csv`: This is the species tree with events labeled. This outputs the Newick representation used in Figure S3 (for both ReconcILS and ete3).
- Timing files are not relevant to our experiments but record timing for ReconcILS and writing output to the files.
- `ZeroCol_ASTRAL_ML_<SCO/all_paralogy>_MIN27_pruned_Gene_trees.csv`: These are the gene trees that produced the particular results. This is important because the input gene trees are non-binary and have 29 species (including outgroups). We have pruned the species, converted them to letters with the key provided above, and resolved the non-binary nodes at random. Hence, this file is the true input to ReconcILS.

### Simulated Dataset
By the end of the run, you will have three types of files:

1. **Reconciliation Results**
   - `<folder_name>_1_results.csv`: Records all the gene trees, with event numbers (i.e., number of duplications, NNI, and loss) required for reconciling each gene tree with the species tree for (True_process <DupCoal>, ReconcILS, DLCpar, and ete3).

2. **Timing**
   - `<folder_name>_1_time_result.csv`: Timing for DLCPar and ReconcILS.

3. **Errors**
   - `<folder_name>_error.txt`: Any error encountered.

## Other Files:
1. **concord.cf.stat**: Gene tree concordance factor for SCO trees.
2. **Plot.r**: Used to plot results in Figure 5 in the main text. Please copy the labeled Newick string from `ZeroCol_ASTRAL_ML_<SCO/all_paralogy>_MIN27_pruned_Gene_trees.csv` to `newick_str` (line 10) and change the `Normalizing_factor` (line 9). SCO: 1820, All_paralogy: 11555.
3. **Produce_results.ipynb**: Jupyter notebook used to produce Figure 4 in the main text.



## MAP Used:
```python
species_to_letters = {
    'Carlitosyrichta': 'A',
    'Cebuscapucinus': 'B',
    'Cercocebusatys': 'C',
    'Macacafascicularis': 'D',
    'Macacanemestrina': 'E',
    'Theropithecusgelada': 'F',
    'Papioanubis': 'G',
    'Macacamulatta': 'H',
    'Mandrillusleucophaeus': 'I',
    'Chlorocebussabaeus': 'J',
    'Colobusangolensis': 'K',
    'Piliocolobustephrosceles': 'L',
    'Rhinopithecusbieti': 'M',
    'Rhinopithecusroxellana': 'N',
    'Gorillagorilla': 'O',
    'Homosapiens': 'P',
    'Panpaniscus': 'Q',
    'Pantroglodytes': 'R',
    'Pongoabelii': 'S',
    'Nomascusleucogenys': 'T',
    'Saimiriboliviensis': 'U',
    'Aotusnancymaae': 'V',
    'Callithrixjacchus': 'W',
}
