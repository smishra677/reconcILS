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

genetree = PhyloTree('(((C,B),C),(((((C,B),C),((((C,B),B),(B,C)),B)),((A,A),((((A,A),(A,A)),A),A))),((((A,A),A),(((((C,B),B),C),(C,B)),B)),((((A,A),A),(((((C,C),B),(C,B)),((C,B),((C,C),((C,B),C)))),B)),((((B,C),B),B),(((C,B),((C,C),B)),(B,C)))))));')
sptree = PhyloTree('(A,(B,C));')

recon_tree_ete, events = genetree.reconcile(sptree)
recon_tree_ete.show()

print(Counter(sp_event_ete3(recon_tree_ete)))