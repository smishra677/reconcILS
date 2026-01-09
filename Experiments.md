# Dataset

## Input

* Clone `reconcILS` from the main branch.
* Download `Experiments_V5` and extract it into `reconcILS/example/Experiments_V5`.
* `Data_with_duplication_and_loss/Data_with_duplication_3_0.03_0.03`
* `Data_with_duplication_and_loss/Dupcoal_output_40_0.03_0.03`
* `Data_with_duplication_and_loss/simphy_output_40_0.0057_0.0047`
* `Data_with_duplication_and_loss/simphy_output_40_0.0057_0.0047`
* `Data_with_duplication_and_loss/simphy_output_40_0.0057_0.0047`
* `Data_with_duplication_and_loss/Data_with_duplication_40_species_5k_0.0003_0.0003_1`
* `Data_with_duplication_and_loss/Data_with_duplication_40_species_2k_0.0003_0.0003`
* `Data_with_duplication_and_loss/Data_with_duplication_fungi_1_0_1_0.015_0.03`

### Folders:

* These folders contain gene trees simulated by DupCoal.
* The species tree used is inside the folder; for example, "species_3_with_branch_lengths.tre".
* The command given to DupCoal is also found inside the folder, named "command.command".

### Files

* **ZeroCol_ASTRAL_ML_ALLPARALOGS_MIN27.tre, ZeroCol_ASTRAL_ML_SCO_MIN27.tre**: These are gene trees from Smith et al., 2020.
* **species_3_with_branch_lengths.tre, species_3_without_branch_lengths.tre, species_40_with_branch_lengths.tre, species_30_without_branch_lengths.tre, species_tree_primates.tre, sp_tree_pruned.tre**: These are the species trees used in our experiments.

  1. **species_3_with_branch_lengths.tre**: Used in simulated experiments (`Data_with_duplication_and_loss/Data_with_duplication_3_0.03_0.03`). The recipe can be found in the DupCoal recipe file.
  2. **species_40_with_branch_lengths.tre**: Used in simulated experiments (`Data_with_duplication_and_loss/Dupcoal_output_40_0.03_0.03`, Figure 4). The recipe can be found in the DupCoal recipe file.
  3. **sp_tree_pruned.tre**: This is the species tree derived from "species_tree_primates.tre" after pruning outgroups and converting species into keys provided in the map.
  4. **fungi_scaled_0.1.tre**: Used in simulated experiments (`Data_with_duplication_and_loss/Data_with_duplication__fungi_1_0_1_0.015_0.03`). The recipe can be found in the DupCoal recipe file.

# Experiments

## Requirements

1. Extract the `Experiments_V4` folder into `./reconcILS`.
2. Install `ete3` and `DLCpar`.

## Single Copy Orthologs / All Paralogy Dataset

3. Run the appropriate script: `python ./run_<SCO/all_paralogy>.py`.

## Simulated Dataset

4. Run the following:

   * `python trigger_simulated_data_3_sp.py`
   * `python trigger_simulated_data_16_sp.py`
   * `python test_reconcILS_default_40.py`
   * `python test_reconcILS_1_1_5_40.py`
   * `python test_reconcILS_3_4_1_40.py`
     *(DLCpar is run for only 3 and 16 species.)*

5. Add noise:

   * `python add_noise.py --file <Path to gene tree folder (all the simulated gene trees for 40 species)> --error '<Average number of NNIs to add to the gene tree>'`

6. Run the following on errored data:

   * `python test_reconcILS_default_40_error_arg.py --file <"Error file">`

# Outputs

### Single Copy Orthologs / All Paralogy Dataset

By the end of the run, the program will output three types of files:

1. **Reconciliation Results**

   * `bio_result_ZeroCol_ASTRAL_ML_<SCO/all_paralogy>_MIN27_pruned.csv`
   * `labeled_L_ZeroCol_ASTRAL_ML_<SCO/all_paralogy>_MIN27_pruned.csv`

2. **Timing**

   * `timing_ZeroCol_ASTRAL_ML_<SCO/all_paralogy>_MIN27_pruned.csv`
   * `timing_ZeroCol_ASTRAL_ML_<SCO/all_paralogy>_MIN27_intro_1_write_pruned.csv`

3. **Output After Pruning Gene Trees**

   * `ZeroCol_ASTRAL_ML_<SCO/all_paralogy>_MIN27_pruned_Gene_trees.csv`

**Output Explanation:**

* `bio_result_ZeroCol_ASTRAL_ML_<SCO/all_paralogy>_MIN27_pruned.csv`: Records all the gene trees, along with event numbers (i.e., the number of duplications, NNI, and losses) required for reconciling each gene tree with the species tree for both reconcILS and ete3.
* `labeled_L_ZeroCol_ASTRAL_ML_SCO_MIN27_pruned.csv`: This is the species tree with events labeled. This file outputs the Newick representation used in Figure S9 (for both reconcILS and ete3).
* Timing files are not relevant to our experiments but record the timing for reconcILS and writing output to the files.
* `ZeroCol_ASTRAL_ML_<SCO/all_paralogy>_MIN27_pruned_Gene_trees.csv`: These are the gene trees that produced the particular results. This is important because the input gene trees are non-binary and include 29 species (including outgroups). We have pruned the species, converted them to letters using the key provided above, and resolved the non-binary nodes randomly. Hence, this file is the true input to reconcILS.

### Simulated Dataset

By the end of the run, you will have four types of files inside the data folder:

1. **Reconciliation Results**

   * `reconcILS_results.csv`: Records all the gene trees, along with event numbers (i.e., the number of duplications, NNI, and losses) required for reconciling each gene tree with the species tree (True_process <DupCoal>, reconcILS, DLCpar, and ete3).

2. **Timing**

   * `reconcILS_time_result.csv`: Timing for ete3 and reconcILS.

3. **Errors**

   * `error.txt`: Any errors encountered.

4. **Memory**

   * `reconcILS_memory_result.csv`: Memory usage for ete3 and reconcILS.

## Other Files:

1. **concord.cf.stat**  
   Gene tree concordance factor for SCO trees.

2. **Plot.r**  
   Script used to plot the results shown in Figure 5 of the main text.  
   Copy the labeled Newick string from  
   `ZeroCol_ASTRAL_ML_<SCO/all_paralogy>_MIN27_pruned_Gene_trees.csv`  
   into `newick_str` (line 10), and set the `Normalizing_factor` as follows:  
   - SCO: `1820`  
   - All_paralogy: `11555`

3. **Produce_results-40_species_V5.ipynb**  
   Jupyter notebook used to generate Figures 4, S6, and S8 in the main text.  
   This notebook also includes:
   - The hemiplasy significance test  
   - Time and memory plots (Figure S4)  
   - Table S1 and Table S2

4. **Produce_results-3_species_V5.ipynb**  
   Jupyter notebook used to produce Figure S7 in the main text.

5. **bin/**  
   Directory containing GCF files for SimPhy and DupCoal gene trees.

6. **produce_tree.r**  
   Script to generate a random tree using the `ape` package.

7. **scale.py**  
   Script used to scale tree branch lengths.

8. **plot_concord_primates.py**  
   Script used to plot and compare correlations for:
   - NNI vs. GCF from Smith et al.
   - All-paralogy NNI vs. SCO NNI
   - ETE3 duplications vs. GCF from Smith et al.  
   Uncomment the relevant sections of the code to generate each comparison.

9. **Data_with_duplication_loss/**  
   Contains gene trees generated by SimPhy and DupCoal, along with all associated results.

10. **Primate_results/**  
    Contains the following files:
    - `all_paralogy.jpg`
    - `all_paralogy.png`
    - `SCO.jpg`
    - `SCO.png`
    - `bio_result_ZeroCol_ASTRAL_ML_ALLPARALOGS_MIN27_pruned.csv`
    - `bio_result_ZeroCol_ASTRAL_ML_SCO_MIN27_pruned_1.csv`
    - `labeled_L_ZeroCol_ASTRAL_ML_ALLPARALOGS_MIN27_pruned.csv`
    - `labeled_L_ZeroCol_ASTRAL_ML_SCO_MIN27_pruned_1.csv`
    - `timing_ZeroCol_ASTRAL_ML_ALLPARALOGS_MIN27_intro_1_write_pruned.csv`
    - `timing_ZeroCol_ASTRAL_ML_SCO_MIN27_intro_1_write_pruned_1.csv`


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

```
