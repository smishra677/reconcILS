from Bio import Phylo
from io import StringIO
from Bio.Phylo.TreeConstruction import *
from ete3 import Tree,TreeNode
from  ete3 import phylo
from ete3 import PhyloTree
from mock_function import get_reconciled_tree,get_reconciled_tree_zmasek
from itertools import permutations
from ete3 import NCBITaxa
import budgitree
from collections import deque
import Bio
#trees = Phylo.read("example.dnd", "newick")

def read_tree():
    trees = Tree("(A, (B, C));")
    return trees

def generator_(species, n):
    for perm in permutations(species, n):
        yield tuple([tuple(perm), tuple(s for s in species if s not in perm)])

def permute(recon_tree,genetree):
    recon_tree_dict={}
    eve_dic={}
    aa = (genetree.get_leaves())
    li=[]
    for j in aa:
        li.append(str(j)[3:])

    

    ge_list =return_perbutation(li)

    for ge in ge_list:
        try:
            genetree1 = PhyloTree(ge)
            
            
            genetree1.resolve_polytomy(recursive=False)
            if len(genetree1.get_leaves())>2:

                recon_tree1, events1 = get_reconciled_tree(genetree1,sptree,[],visited)
                recon_tree_dict[ge]= recon_tree1
                eve_dic[ge]=events1
        except:
            continue
    return recon_tree_dict,eve_dic


def eve_print(events):
    for ev in events:
        if ev.etype == "S":
            print('ORTHOLOGY RELATIONSHIP:', ','.join(ev.inparalogs), "<====>", ','.join(ev.orthologs))
        elif ev.etype == "L":
            print('Loss RELATIONSHIP:', ','.join(ev.inparalogs), "<====>", ','.join(ev.outparalogs))
        else:
            print('PARALOGY RELATIONSHIP:', ','.join(ev.inparalogs), "<====>", ','.join(ev.outparalogs))
            '''
            print('OR')
            aru= [ '((A,C),B);', '((A,B),C);', '((B,C),A);']
            for  i in aru:
                genetree1 = PhyloTree(i)
                recon_tree1, events1 = get_reconciled_tree(genetree1,sptree,[],visited)
                for ev1 in events1:
                    if ev1.etype == "S":
                        print('ORTHOLOGY RELATIONSHIP:', ','.join(ev1.inparalogs), "<====>", ','.join(ev1.orthologs))
                    elif ev.etype == "L":
                        print('Loss RELATIONSHIP:', ','.join(ev1.inparalogs), "<====>", ','.join(ev1.outparalogs))
                    else:
                        print('PARALOGY RELATIONSHIP:', ','.join(ev1.inparalogs), "<====>", ','.join(ev1.outparalogs))
            '''


def return_perbutation(species):
    newicks= []
    for n in range(1, len(species)):
        for p in generator_(species, n):
            newicks.append('((' + ', '.join('{})'.format(', '.join(k)) for k in p)+';'.format(''))

    return newicks

gene_tree_nw ='(((A,D),C),B);'
species_tree_nw = '(((A,B),C),D);'

tr= '(((A,D),C),B);'
sp ='(((A,B),C),D);'
species = ["A","B","C"]




visited={}

genetree = PhyloTree(gene_tree_nw)
sptree = PhyloTree(species_tree_nw)
recon_tree, events = genetree.reconcile(sptree)

eve_print(events)
'''

#print(genetree.robinson_foulds(sptree))


cost ={'D':2,'NNI':1}
recon_tree, events = get_reconciled_tree(genetree,sptree,[],visited)
eve_print(events)

recon_tree.show()

exit()
print(cost)


print(recon_tree.evoltype)
recon_tree_dict ={}
eve_dic= {}
if recon_tree.evoltype=='D':
    recon_tree_dict,eve_dic=permute(recon_tree,genetree)

eve={}

if len(recon_tree_dict)>1:
    for r in recon_tree_dict:
        eve_l=[]
        for node in recon_tree_dict[r].traverse(strategy="postorder"):
            for i in node.children:
                if len(i.get_leaves())>1:
                    eve_l.append(i.evoltype)
                        

        eve[r]=eve_l


for i in eve:
    print(i)
    print(species_tree_nw )
    eve_print(eve_dic[i])
    print(eve[i])
    #recon_tree_dict[i].show()
    print('###############################################################')

#genetree.show()

#recon_tree.show()
'''