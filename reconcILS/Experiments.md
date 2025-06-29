# Dataset

## Input 
- Clone `reconcILS` from the main branch.
- Download `Experiments_V4` and extract it into `reconcILS/Experiments_V4`.
- Data_40_with_duplication_loss/Data_with_duplication_3_0.03_0.03
- Data_40_with_duplication_loss/Dupcoal_output_40_0.03_0.03
- Data_40_with_duplication_loss/simphy_output_40_0.0057_0.0047

### Folders: 
- These folders contain gene trees simulated by dupcoal.
- The species tree used is inside the folder, for example
 "species_3_with_branch_lengths.tre".
- The command given to DupCoal is also found inside the folder, named "command.command".

### Files
- **ZeroCol_ASTRAL_ML_ALLPARALOGS_MIN27.tre, ZeroCol_ASTRAL_ML_SCO_MIN27.tre**: These are gene trees from Smith et al., 2020.
- **species_3_with_branch_lengths.tre, species_3_without_branch_lengths.tre, species_40_with_branch_lengths.tre, species_30_without_branch_lengths.tre, species_tree_primates.tre, sp_tree_pruned.tre**: These are the species trees used in our experiments.
  1. **species_3_with_branch_lengths.tre**: Used in simulated experiments (Data_40_with_duplication_loss/Data_with_duplication_3_0.03_0.03). You can find the recipe in the DupCoal recipe file.
  2. **species_40_with_branch_lengths.tre**: Used in simulated experiments (Data_40_with_duplication_loss/Dupcoal_output_40_0.03_0.03 (Figure 4)). You can find the recipe in the DupCoal recipe file.
  3. **sp_tree_pruned.tre**: This is the species tree from "species_tree_primates.tre" after pruning outgroups and converting species into keys provided in the map.

# Experiments

## Requirements
1. Extract the Experiments_V4 folder into `./reconcILS`.
2. Install `ete3` and `DLCpar`.

## Single Copy Orthologs / All Paralogy Dataset
3. Run the appropriate script: `python ./run_<SCO/all_paralogy>.py`.

## Simulated Dataset
4. Run the following:
   - `python trigger_simulated_data_3_sp.py`
   - `python test_reconcILS_default_40.py`
   - `python test_reconcILS_1_1_5_40.py`
   - `python test_reconcILS_3_4_1_40.py`
   (DLCpar is run for only 3 species.)

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
- `bio_result_ZeroCol_ASTRAL_ML_<SCO/all_paralogy>_MIN27_pruned.csv`: Records all the gene trees, along with event numbers (i.e., number of duplications, NNI, and losses) required for reconciling each gene tree with the species tree for both reconcILS and ete3.
- `labeled_L_ZeroCol_ASTRAL_ML_SCO_MIN27_pruned.csv`: This is the species tree with events labeled. This outputs the Newick representation used in Figure S9 (for both reconcILS and ete3).
- Timing files are not relevant to our experiments but record the timing for reconcILS and writing output to the files.
- `ZeroCol_ASTRAL_ML_<SCO/all_paralogy>_MIN27_pruned_Gene_trees.csv`: These are the gene trees that produced the particular results. This is important because the input gene trees are non-binary and include 29 species (including outgroups). We have pruned the species, converted them to letters using the key provided above, and resolved the non-binary nodes at random. Hence, this file is the true input to reconcILS.

### Simulated Dataset
By the end of the run, you will have four types of files inside the data folder:

1. **Reconciliation Results**
   - `reconcILS_results.csv`: Records all the gene trees, along with event numbers (i.e., number of duplications, NNI, and losses) required for reconciling each gene tree with the species tree (True_process <DupCoal>, reconcILS, DLCpar, and ete3).

2. **Timing**
   - `reconcILS_time_result.csv`: Timing for ete3 and reconcILS.

3. **Errors**
   - `error.txt`: Any errors encountered.

4. **Memory**
   - `reconcILS_memory_result.csv`: Memory usage for ete3 and reconcILS.

## Other Files:
1. **concord.cf.stat**: Gene tree concordance factor for SCO trees.
2. **Plot.r**: Used to plot results in Figure 5 of the main text. Please copy the labeled Newick string from `ZeroCol_ASTRAL_ML_<SCO/all_paralogy>_MIN27_pruned_Gene_trees.csv` into `newick_str` (line 10) and adjust the `Normalizing_factor` as follows: SCO: 1820, All_paralogy: 11555.
3. **Produce_results-40_species_V4.ipynb**: Jupyter notebook used to produce Figure 4, S6, S8 in the main text. You can also find the significance test for hemiplasy in this notebook, as well as Time and Memory plots (Figure S4).
4. **Produce_results-3_species_V4.ipynb**: Jupyter notebook used to produce Figure S7 in the main text.
5. **bin**: Folder containing gcf files for SimPhy and DupCoal gene trees.
6. **produce_tree.r**: Script to generate a random tree using `ape`.
7. **scale.py**: Script to scale the tree length.
8. **plot_concord_primates.py**: Used to plot and compare correlations for (NNI vs gcf from Smith et al.), (all paralogs NNI vs SCO NNI), and (ETE3 duplications vs gcf from Smith et al.). Uncomment different sections of the code to obtain these.
9. **Data_with_duplication_loss**: Contains folders for gene trees produced by SimPhy and DupCoal along with all results.
10. **Primate_results**: Contains the results for the primate dataset.

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
