from Bio import Phylo
from io import StringIO
from Bio.Phylo.TreeConstruction import *

#trees = Phylo.read("example.dnd", "newick")

def read_tree():
    trees = Phylo.read(StringIO("(A, (B, C), (D, E))"), "newick")
    return trees
scorer = ParsimonyScorer()
searcher = NNITreeSearcher(scorer)._get_neighbors(trees)
print(searcher)
#constructor = ParsimonyTreeConstructor(searcher, trees)
#pars_tree = constructor.build_tree()
#print(pars_tree)