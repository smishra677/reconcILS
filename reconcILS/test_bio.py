import sys
sys.path.append("../")
import utils.Tree as Tree
import copy
import utils.Tally as Tally
import argparse
import utils.ILS as ILS
import utils.readWrite as readWrite
import pickle
import reconcILS

def print_(sp):
    if sp:
        if (sp.taxa==None):
            print(sp.taxa)
        print_(sp.leftChild)
        print_(sp.rightChild)


er_sp= readWrite.readWrite().parse_bio(reconcILS.reconcils().read_trees('./sp.tre'))
#print(reconcILS.reconcils().read_trees('./sp.tre'))
#print(er_sp.to_newick())
#print(er.to_newick())
er= readWrite.readWrite().parse_bio(reconcILS.reconcils().read_trees('./gt.tre'))
from ete3 import PhyloTree


gene_tree_nw =  reconcILS.reconcils().read_trees('./gt.tre')
species_tree_nw = reconcILS.reconcils().read_trees('./sp.tre')
genetree = PhyloTree(gene_tree_nw)
sptree = PhyloTree(species_tree_nw)
recon_tree, events = genetree.reconcile(sptree)

