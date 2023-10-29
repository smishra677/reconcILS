from ete3 import TreeNode
from  ete3 import phylo
from ete3 import PhyloTree
from collections import Counter
from ete3 import NCBITaxa
import os

import sys


def sp_event_ete3(root):
        li=[]
        for node in root.traverse(strategy="postorder"):
            if len(node.children) != 0:
            ##print(dir(root))

                print(li)
                
                li.append(node.evoltype)
            ##print(root.isLeaf),
        return li

genetree = PhyloTree('((CCC_0,BBB_0),CCC_1);')
sptree = PhyloTree("(BBB,CCC);")


recon_tree_ete, events = genetree.reconcile(sptree)
recon_tree_ete.show()

print(Counter(sp_event_ete3(recon_tree_ete)))