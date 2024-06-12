# Example Use

## Species Tree String
Please change the species string on line 18 in `example.py`. Alternatively, you can read a species tree in Newick format from a file as a string.
sp_string='(A,(B,C));'


## Change the Costs of Each Event
These are set with default costs on lines 20-31 in `example.py`:
Duplication_cost = 1.1
ILS_cost = 1.0
Loss_cost = 1.0


## Gene Tree as a List of Strings
Alternatively, you can write a file reader to read it as a string or read the gene tree from an individual file by changing line 55:
list_gene_trees = ['((A,(B,C)),(B,(B,C)))']
