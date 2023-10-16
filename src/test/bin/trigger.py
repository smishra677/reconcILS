import os

import sys
sys.path.append("../../")
from reconcILS import *
from utils import *
import pandas as pd     

from ete3 import TreeNode
from  ete3 import phylo
from ete3 import PhyloTree

from ete3 import NCBITaxa
import os





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

def read_trees(i,folder):
        gene_tre= open(folder+'/rep_'+str(i)+'.tre')
        tr =gene_tre.read().strip().split('\n')
        gene_tre.close()
        return str(tr[0])

def write_trees(i,folder,tree):
        with open(folder+'/rep_1_'+str(i)+'.tre', 'w') as f:
            f.write(tree)
        

import timeit

def dlcpar(gene_folder,i):
    command1= 'python2 ./dlcpar dp -D 1.1 -L 1.0 -C 1.0 -K 1.0 -s ./sp_tree.tre -S ./001/001.mapsl  ./'+gene_folder+'/rep_1_'+str(i)+'.tre'
    command2= 'python2 ./dlcpar events --lct -s ./sp_tree.tre -S ./001/001.mapsl ./'+gene_folder+'/rep_1_'+str(i)+'.tre.dlcdp.lct.tree > result_data'
    os.system(command1)
    os.system(command2)

#sp_string='(((A,(B,C)),D),(E,F));'
sp_string='(A,(B,C));'



gene_folder='10_16'
#gene_tre= open('./output_gene/gene_tree.txt')
#trees =gene_tre.read().strip().split('\n')
#gene_tre.close()


value=0
previous_dict={}
erro=0
fil= open(gene_folder+'_error.txt','w+')
#sp ='(((A,B),C),D);'
dlcparTime=[]
reconcILSTime_iterative=[]
reconcILSTime_recurssive=[]
ete_time=[]
li_gt=[]
lis_rep=[]
for i in range(1000):
    red= readWrite.readWrite()
    sp=red.parse(sp_string)
    if erro==1:
        dic=previous_dict

    previous_dict=dic

    DATA_PATH = "./"
    file_path = os.path.join(DATA_PATH)

   
    try:
        tree= str(read_trees(i,os.path.join(file_path, gene_folder)))
        tr=red.parse(tree)
        print(tr.to_newick())
        write_trees(i,os.path.join(file_path, gene_folder),tr.to_newick())

        try:
            starttime = timeit.default_timer()
            command1= 'python2 ./dlcpar dp -D 1.1 -L 1.0 -C 1.0 -K 1.0 -s ./sp_tree.tre -S ./001/001.mapsl  ./'+gene_folder+'/rep_1_'+str(i)+'.tre'
            command2= 'python2 ./dlcpar events --lct -s ./sp_tree.tre -S ./001/001.mapsl ./'+gene_folder+'/rep_1_'+str(i)+'.tre.dlcdp.lct.tree > result_data'
            os.system(command1)
            os.system(command2)

            dlcparTime.append(timeit.default_timer()-starttime)
            df = pd.read_csv('result_data',header=0, delimiter='\t')

            dic= red.read_log('True Process',i,dic,gene_folder)
            lis_rep.append(i)
            li_gt.append(tr.to_newick())

            dic['Gene_tree']+=[tr.to_newick()]
            dic['Species_Tree']+=[sp_string]
            
      

            
            dic['Process'].append('DLCpar')
            dic['Replicate'].append(str(i))
            dic['Gene_tree'].append(tr.to_newick())
            dic['Species_Tree'].append(sp_string)
            dic=convert(dic,df)

        
          

            
            genetree = PhyloTree(tr.to_newick())
            
            sptree = PhyloTree(sp_string)
  
            
            starttime = timeit.default_timer()
            recon_tree_ete, events = genetree.reconcile(sptree)
            ete_time.append(timeit.default_timer()-starttime)
            dic['Gene_tree']+=[tr.to_newick()]
            dic['Species_Tree']+=[sp_string]
        
            dic= red.Create_pd_ete3('ETE3',i,dict(Counter(red.sp_event_ete3(recon_tree_ete))),dic)


            
        except:
            erro=1
            dic=previous_dict
            print(i)
            dlcparTime=dlcparTime[:-1]
            #li_gt=li_gt[:-1]
            #lis_rep=lis_rep[:-1]
            fil.write("dlcpar "+tree+"   "+str(i)+'\n')
            print(tr.to_newick())
            continue
 
        #print(dic)

       

        try:
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

            starttime1 = timeit.default_timer()
            #reconcILS(tr,sp,sp_copy,sp,[])
            li=iterative_reconcILS(tr,sp,sp_copy,sp,[])

            reconcILSTime_iterative.append(timeit.default_timer()-starttime1)
   
      
            dic['Gene_tree']+=[tr.to_newick()]
            dic['Species_Tree']+=[sp_string]
            #li =red.sp_event(sp,[])

            dic= red.Create_pd('Our_algorithm_ite',i,li,dic)

        except:
            erro=1
            dic=previous_dict


            fil.write("reconcILS_iterative"+"   "+str(i)+"  "+tree+'\n')

        try:
            sp=red.parse(sp_string)
            sp_copy= copy.deepcopy(sp)
            sp_copy.reset()
            sp.reset()
            tr.reset()
            tr.order_gene(sp)

            tr.label_internal()
            sp.label_internal()



            tr.map_gene(sp)
            setCost(sp)

            sp.isRoot=True
            tr.isRoot=True
            sp_copy.isRoot=True

            starttime1 = timeit.default_timer()
            reconcILS(tr,sp,sp_copy,sp,[])

            #li=iterative_reconcILS(tr,sp,sp_copy,sp,[])
            reconcILSTime_recurssive.append(timeit.default_timer()-starttime1)
            
      
            dic['Gene_tree']+=[tr.to_newick()]
            dic['Species_Tree']+=[sp_string]
            li =red.sp_event(sp,[])
            #print(Counter(li))

            dic= red.Create_pd('Our_algorithm_recu',i,li,dic)

        except:
            erro=1
            dic=previous_dict


            fil.write("reconcILS"+"   "+str(i)+"  "+tree+'\n')


    


    except:
        erro=1
        dic=previous_dict
        print(i)
        #dlcparTime=dlcparTime[:-1]
        fil.write("Empty: "+str(i)+'\n')
        print(tr.to_newick())
        continue

fil.close()

            
if erro==1:
    dic=previous_dict
for i in dic:
    print(len(dic[i]))
    
print(dic)
df = pd.DataFrame(dic)
print(df)
print(len(lis_rep),len(li_gt))
time_dic={'Replicate':lis_rep,'Gene_tree':li_gt,'reconcILSTime_iterative':reconcILSTime_iterative,'reconcILSTime_recurssive':reconcILSTime_recurssive,'dlcparTime':dlcparTime,'ete_time':ete_time}
print(time_dic)
print(len(reconcILSTime_iterative))
print(len(reconcILSTime_recurssive))

print(len(dlcparTime))
df1 = pd.DataFrame(time_dic)



df.to_csv(gene_folder+'_result.csv', index=False)

df1.to_csv(gene_folder+'_time_result.csv', index=False)


        


        


   
