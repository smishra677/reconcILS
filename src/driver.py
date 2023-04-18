from Bio import Phylo
from io import StringIO
from Bio.Phylo.TreeConstruction import *
from ete3 import Tree
from  ete3 import phylo
from ete3 import PhyloTree
from mock_function import get_reconciled_tree,get_reconciled_tree_zmasek
#trees = Phylo.read("example.dnd", "newick")

def read_tree():
    trees = Tree("(A, (B, C));")
    return trees


gene_tree_nw = '((A_001,B_001),(A_002,B_002));'
species_tree_nw = "(A_0,B_0);"
genetree = PhyloTree(gene_tree_nw)
sptree = PhyloTree(species_tree_nw)

#recon_tree, events = get_reconciled_tree_zmasek(genetree,sptree)
recon_tree = get_reconciled_tree_zmasek(genetree,sptree)



#phylo.reconciliation. get_reconciled_tree(t,t,[])
#scorer = ParsimonyScorer()
#searcher = NNITreeSearcher(scorer)._get_neighbors(trees)
#
# print(searcher)
#constructor = ParsimonyTreeConstructor(searcher, trees)
#pars_tree = constructor.build_tree()
#print(pars_tree)

'''

for ev in events:
    if ev.etype == "S":
        print('ORTHOLOGY RELATIONSHIP:', ','.join(ev.inparalogs), "<====>", ','.join(ev.orthologs))
    elif ev.etype == "D":
        print('PARALOGY RELATIONSHIP:', ','.join(ev.inparalogs), "<====>", ','.join(ev.outparalogs))
    else:
        print('Loss RELATIONSHIP:', ','.join(ev.inparalogs), "<====>", ','.join(ev.outparalogs))
'''
print(recon_tree.evoltype)
for node in recon_tree.traverse(strategy="postorder"):
    for i in node.children:
        if len(i.get_leaves())>1:
            print(i.evoltype)
#genetree.show()

#recon_tree.show()