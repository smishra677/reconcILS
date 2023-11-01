import os
import copy

import sys
sys.path.append("../")
from reconcILS import *
from utils import *
import pandas as pd 

def read_trees(i,folder):
        gene_tre= open(folder+'/rep_'+str(i)+'.tre')
        tr =gene_tre.read().strip().split('\n')
        gene_tre.close()
        return str(tr[0])


sp_string='(A,(B,C));' # species Tree string
dic={'Process':[],'Replicate':[],'Gene_tree':[],'Species_Tree':[],'Duplication':[],'NNI':[],'Loss':[]}

DATA_PATH = "./" # path to the gene tree folder
file_path = os.path.join(DATA_PATH)

gene_folder='10_30'

n =1000 # number of gene tree replicate in the folder
for i in range(1):
    reco =reconcils()
    red= readWrite.readWrite()

    tree= str(read_trees(i,os.path.join(file_path, gene_folder)))
    tr=red.parse(tree)

    sp=red.parse(sp_string)
    sp_copy= copy.deepcopy(sp)
    sp_copy.reset()
    sp.reset()
    

    tr.order_gene(sp)
    tr.label_internal()
    sp.label_internal()
    tr.map_gene(sp)
    reco.setCost(sp)
    sp.isRoot=True
    tr.isRoot=True
    sp_copy.isRoot=True


    reco.gene_tree= copy.deepcopy(tr)
   
    reco.reconcILS(tr,sp,sp_copy,sp,[])

            
      
    dic['Gene_tree']+=[red.to_newick(reco.gene_tree)]
    dic['Species_Tree']+=[sp_string]
    li =red.sp_event(sp,[])
    print(li)

    dic= red.Create_pd('reconcILS',i,li,dic)

df = pd.DataFrame(dic)
print(df)


df.to_csv(gene_folder+'_result.csv', index=False)


