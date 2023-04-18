from Bio import Phylo
from io import StringIO
from Bio.Phylo.TreeConstruction import *
from ete3 import Tree,TreeNode
from  ete3 import phylo
from ete3 import PhyloTree
from mock_function import get_reconciled_tree,get_reconciled_tree_zmasek
from itertools import permutations


import Bio
#trees = Phylo.read("example.dnd", "newick")

def read_tree():
    trees = Tree("(A, (B, C));")
    return trees

def generator_(species, n):
    for perm in permutations(species, n):
        yield tuple([tuple(perm), tuple(s for s in species if s not in perm)])





def return_perbutation(species):
    newicks= []
    for n in range(1, len(species)):
        for p in generator_(species, n):
            newicks.append('((' + ', '.join('{})'.format(', '.join(k)) for k in p)+');'.format(''))
    
    return newicks

gene_tree_nw = '((A,C),B);'
species_tree_nw = "(A, (B, C));"

species = ["A","B","C"]

print(return_perbutation(species))
visited={}

genetree = PhyloTree(gene_tree_nw)
sptree = PhyloTree(species_tree_nw)
cost ={'D':2,'NNI':1}
recon_tree, events,cost = get_reconciled_tree(genetree,sptree,[],cost,visited)

print(cost)
#print(Tree(recon_tree. get_speciation_trees()[2]))
#recon_tree = get_reconciled_tree_zmasek(genetree,sptree)
scorer = ParsimonyScorer()

searcher = NNITreeSearcher(scorer)
constructor = ParsimonyTreeConstructor(searcher, gene_tree_nw)
#for i in constructor:
#print(i)
#print(searcher)

#phylo.reconciliation. get_reconciled_tree(t,t,[])
#scorer = ParsimonyScorer()
#searcher = NNITreeSearcher(scorer)._get_neighbors(trees)
#
# print(searcher)
#constructor = ParsimonyTreeConstructor(searcher, trees)
#pars_tree = constructor.build_tree()
#print(pars_tree)



for ev in events:
    if ev.etype == "S":
        print('ORTHOLOGY RELATIONSHIP:', ','.join(ev.inparalogs), "<====>", ','.join(ev.orthologs))
    elif ev.etype == "D":
        print('PARALOGY RELATIONSHIP:', ','.join(ev.inparalogs), "<====>", ','.join(ev.outparalogs))
    else:
        print('Loss RELATIONSHIP:', ','.join(ev.inparalogs), "<====>", ','.join(ev.outparalogs))
print(recon_tree.evoltype)
for node in recon_tree.traverse(strategy="postorder"):
    for i in node.children:
        if len(i.get_leaves())>1:
            print(i.evoltype)
#genetree.show()

recon_tree.show()