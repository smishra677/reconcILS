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
      
        print(sp.get_species())
        for k in sp.children:
             print_(k)



from ete3 import PhyloTree

# Loads a gene tree and its corresponding species tree. Note that
# species names in sptree are the 3 firs letters of leaf nodes in
# genetree.
gene_tree_nw = '(A,((B,((C,((D,E),F)),((G,((((((H,((I,J),K)),(L,M)),N),O),((P,Q),(R,S))),(((T,(U,(V,W))),X),Y))),(Z,(AA,BB))))),CC));'
species_tree_nw = '((((((((U,(W,V)),T),X),Y),(((((N,H),(K,L)),((M,I),J)),O),((S,R),(Q,P)))),((Z,G),(BB,AA))),(((E,D),F),((B,CC),A))),C);'
genetree = PhyloTree(gene_tree_nw)
sptree = PhyloTree(species_tree_nw)


#                    /-Dme_001
#          /--------|
#         |          \-Dme_002
#         |
#         |                              /-Cfa_001
#         |                    /--------|
#---------|                   |          \-Mms_001
#         |          /--------|
#         |         |         |                    /-Hsa_001
#         |         |         |          /--------|
#         |         |          \--------|          \-Ptr_001
#          \--------|                   |
#                   |                    \-Mmu_001
#                   |
#                   |          /-Ptr_002
#                    \--------|
#                             |          /-Hsa_002
#                              \--------|
#                                        \-Mmu_002
#
# Let's reconcile our genetree with the species tree
recon_tree, events = genetree.reconcile(sptree)
# a new "reconcilied tree" is returned. As well as