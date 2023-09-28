import os

import sys
#sys.path.append("../")
from reconcILS import *


import pandas as pd     



def read_trees(i,folder):
     gene_tre= open('./'+folder+'/rep_'+str(i)+'.tre')
     tr =gene_tre.read().strip().split('\n')
     gene_tre.close()
     return str(tr[0])

def write_tree(i,folder,tree):
    with open('./'+folder+'/rep_1_'+str(i)+'.tre', 'w') as f:
        f.write(tree)



dic={'Process':[],'Replicate':[],'Gene_tree':[],'Species_Tree':[],'Duplication':[],'NNI':[],'DLCILS':[],'Loss':[],'Hemiplasy':[],'RHemiplasy':[]}
    
from collections import Counter

def convert(df,df1):
	df['Loss'].append(sum(list(df1['loss'])))
	df['DLCILS'].append(sum(list(df1['coal'])))
	df['Duplication'].append(sum(list(df1['dup'])))
	df['Hemiplasy'].append(0)
	df['NNI'].append(0)
	df['RHemiplasy'].append(0)
	
	return df




sp_string='(A,(B,C))'



gene_folder='9_19_1000'
#gene_tre= open('./output_gene/gene_tree.txt')
#trees =gene_tre.read().strip().split('\n')
#gene_tre.close()


value=0
previous_dict={}
erro=0
fil= open(gene_folder+'_error.txt','w+')
#sp ='(((A,B),C),D);'
for i in range(1000):
    sp=parse(sp_string)
    if erro==1:
        dic=previous_dict

    previous_dict=dic
    try:
        tree= str(read_trees(i,gene_folder))
        tr=parse(tree)
        write_tree(i,gene_folder,to_newick(tr))

        command1= 'python2 ./dlcpar dp -D 1.1 -L 1.0 -C 1.0 -K 1.0 -s ./sp_tree_large.tre -S ./001/001.mapsl  ./'+gene_folder+'/rep_1_'+str(i)+'.tre'
        command2= 'python2 ./dlcpar events --lct -s ./sp_tree_large.tre -S ./001/001.mapsl ./'+gene_folder+'/rep_1_'+str(i)+'.tre.dlcdp.lct.tree > result_data'
        os.system(command1)
        os.system(command2)
        df = pd.read_csv('result_data',header=0, delimiter='\t')
        dic= read_log('True Process',i,dic,gene_folder)
        
        print(dic)
        dic['Gene_tree']+=[to_newick(tr)]
        dic['Species_Tree']+=[sp_string]
        
  

        
        dic['Process'].append('DLCpar')
        dic['Replicate'].append(str(i))
        dic['Gene_tree'].append(to_newick(tr))
        dic['Species_Tree'].append(sp_string)
        dic=convert(dic,df)
        print(dic)

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
        try:
            reconcILS(tr,sp,sp_copy,sp)
            
      
            dic['Gene_tree']+=[to_newick(tr)]
            dic['Species_Tree']+=[sp_string]
            li =sp_event(sp,[])
            #print(Counter(li))
            print(li)
            dic= Create_pd('Our_algorithm',i,li,dic)
            print(dic)
        except:
            erro=1
            dic=previous_dict
            fil.write("reconcILS"+"   "+str(i)+"  "+tree+'\n')





    except:
        erro=1
        dic=previous_dict
        print(i)
        fil.write(tree+"   "+str(i)+'\n')
        print(to_newick(tr))
        continue

fil.close()

                                        
	
            
if erro==1:
    dic=previous_dict
for i in dic:
    print(len(dic[i]))
    
print(dic)
df = pd.DataFrame(dic)
print(df)



df.to_csv(gene_folder+'_result.csv', index=False)


        


   
