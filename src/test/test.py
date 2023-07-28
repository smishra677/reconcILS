from ete3 import TreeNode
from  ete3 import phylo
from ete3 import PhyloTree

from ete3 import NCBITaxa
import os

import sys
sys.path.append("../")
from reconcILS import *



import pandas as pd     




dic={'Process':[],'Replicate':[],'Gene_tree':[],'Species_Tree':[],'Duplication':[],'NNI':[],'Loss':[]}
    
from collections import Counter


sp_string='(A,(B,C));'
'''
gene_tree='((A,C),B);'
tr=parse(gene_tree)

sp=parse(sp_string)
sp_copy= copy.deepcopy(sp)
sp_copy.reset()

tr.order_gene(sp)

tr.label_internal()
sp.label_internal()



tr.map_gene(sp)

setCost(sp)
sp.isRoot=True
tr.isRoot=True
sp_copy.isRoot=True
reconcILS(tr,sp,sp_copy,sp)
#print('######################33')
#print(to_newick(sp))

li =sp_event(sp,[])

print(Counter(li))
exit()
'''

gene_folder='7_25'
#gene_tre= open('./output_gene/gene_tree.txt')
#trees =gene_tre.read().strip().split('\n')
#gene_tre.close()


value=0
previous_dict={}
erro=0
#sp ='(((A,B),C),D);'
for i in range(100):


    sp=parse(sp_string)
    if erro==1:
        dic=previous_dict

    previous_dict=dic

    try:
        tree= str(read_trees(i,gene_folder))
        #print('tree',tree)
        tr=parse(tree)
        tr=parse(to_newick(tr))


        


        sp_copy= copy.deepcopy(sp)
        sp_copy.reset()

        tr.order_gene(sp)

        tr.label_internal()
        sp.label_internal()



        tr.map_gene(sp)
        setCost(sp)

        sp.isRoot=True
        tr.isRoot=True
        sp_copy.isRoot=True
        reconcILS(tr,sp,sp_copy,sp)

        dic['Gene_tree']+=[to_newick(tr)]
        dic['Species_Tree']+=[sp_string]
        
        dic['Gene_tree']+=[to_newick(tr)]
        dic['Species_Tree']+=[sp_string]    

        dic= read_log('True Process',i,dic,gene_folder)

    
        
        #print('######################33')
        #print(to_newick(sp))

        li =sp_event(sp,[])
        #print(Counter(li))

        dic= Create_pd('Our_algorithm',i,dict(Counter(li)),dic)

        dic['Gene_tree']+=[to_newick(tr)]
        dic['Species_Tree']+=[sp_string]


        genetree = PhyloTree(to_newick(tr))
        sptree = PhyloTree(sp_string)

        recon_tree_ete, events = genetree.reconcile(sptree)

        dic= Create_pd_ete3('ETE3',i,dict(Counter(sp_event_ete3(recon_tree_ete))),dic)
    except:
        erro=1
        print(to_newick(tr))

        continue
#print(dic)
if erro==1:
    dic=previous_dict
for i in dic:
    print(len(dic[i]))
df = pd.DataFrame(dic)

df.to_csv('result_7_14_1.csv', index=False)
'''

Initial_multiple_mapping=1
#sp.find_cost(tr,0)
new_topo,cost =(ILS(tr,sp,sp_copy,Initial_multiple_mapping))
#print(to_newick(new_topo))
#print(Initial_multiple_mapping-cost)




#tr.map_recon(recon)
#print(to_newick(recon))
new_tr= copy.deepcopy(tr)
new_tr.reset()
#print(new_tr.optimize_cost(recon,recon))
#recon.tag_loss()

#print(to_newick(recon))

'''