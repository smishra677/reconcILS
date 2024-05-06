import os
import copy
from ete3 import PhyloTree  
from collections import Counter
import sys
sys.path.append("../")
from reconcILS import *
from utils import *
import pandas as pd 
import time
from concurrent.futures import ThreadPoolExecutor, TimeoutError
def read_trees(i,folder):
        gene_tre= open(folder+'/rep_'+str(i)+'.tre')
        tr =gene_tre.read().strip().split('\n')
        gene_tre.close()
        return str(tr[0])

def get_edges(sp):
        return_it=[(repr(list(sorted(sp.taxa))),' to ',repr(list(sorted(sp.taxa))))]
        stack = [sp]
        while stack:
            curr=stack.pop()

            if curr.leftChild:
                return_it.append((repr(list(sorted(curr.taxa))),' to ',repr(list(sorted(curr.leftChild.taxa)))))
                stack.append(curr.leftChild)
            if curr.rightChild:
                return_it.append((repr(list(sorted(curr.taxa))),' to ',repr(list(sorted(curr.rightChild.taxa)))))
                stack.append(curr.rightChild)
        return return_it         

def edge_to_event(sp,dic):
        stack = [sp]
        sp.sp_ev_list+=[[dic[(repr(list(sorted(sp.taxa))),' to ',repr(list(sorted(sp.taxa))))],'Up']]
        while stack:
            curr=stack.pop()
            if curr.parent:
                curr.sp_ev_list+= [[dic[(repr(list(sorted(curr.parent.taxa))),' to ',repr(list(sorted(curr.taxa))))],'Up']]
            if curr.leftChild:
                stack.append(curr.leftChild)
            if curr.rightChild:
                stack.append(curr.rightChild)


def generate_species_codes(mapping,tree):
    for species, letter in mapping.items():
        tree = tree.replace(species, letter)
    return tree
species_to_letters = {
    'Galeopterusvariegatus': 'A',
    'Tupaiachinensis': 'B',
    'Carlitosyrichta': 'C',
    'Microcebusmurinus': 'D',
    'Propithecuscoquereli': 'E',
    'Otolemurgarnettii': 'F',
    'Cebuscapucinus': 'G',
    'Cercocebusatys': 'H',
    'Macacafascicularis': 'I',
    'Macacanemestrina': 'J',
    'Theropithecusgelada': 'K',
    'Papioanubis': 'L',
    'Macacamulatta': 'M',
    'Mandrillusleucophaeus': 'N',
    'Chlorocebussabaeus': 'O',
    'Colobusangolensis': 'P',
    'Piliocolobustephrosceles': 'Q',
    'Rhinopithecusbieti': 'R',
    'Rhinopithecusroxellana': 'S',
    'Gorillagorilla': 'T',
    'Homosapiens': 'U',
    'Panpaniscus': 'V',
    'Pantroglodytes': 'W',
    'Pongoabelii': 'X',
    'Nomascusleucogenys': 'Y',
    'Saimiriboliviensis': 'Z',
    'Aotusnancymaae': 'AA',
    'Callithrixjacchus': 'AB',
    'Musmusculus': 'AC'
}

sp_string='(D,(C,(A,B)));'

dic={'Process':[],'Replicate':[],'Gene_tree':[],'Species_Tree':[],'Duplication':[],'NNI':[],'Loss':[]}
DATA_PATH = "./no_dup" # path to the gene tree folder
file_path = os.path.join(DATA_PATH)
time_it=[]

time_it2=[]

reco =reconcils()
red= readWrite.readWrite()
sp_large=red.parse(sp_string)

sp_str_ref=sp_large.to_newick()

#sp_str_ref = generate_species_codes(species_to_letters,sp_str_ref)


sp_large=sp_large.parse(sp_str_ref)
sp_large.label_internal()
red.write_introgression(sp_large)


import gc


for ili in range(1000):
                print('##########################################################################################################################')
                print(ili)
                print('--------------------------------------------------------------------------------------------------------------------------')

                reco =reconcils()
                red= readWrite.readWrite()
        
                op= open(file_path+'/rep_'+str(ili)+'.tre')
        
                tree= str(op.read())

                if len(tree)==0:
                        continue

                #tree = tree.replace('e-', '0')

                
                
                tr= red.parse(red.to_newick(red.parse(tree)))

                sp=red.parse(sp_string)
                #dic= red.read_log('True Process',i,dic,gene_folder)
                tr_str_ref=tr.to_newick()
                sp_str_ref=sp.to_newick()
                #tr_str_ref = generate_species_codes(species_to_letters,tr_str_ref)
                #sp_str_ref = generate_species_codes(species_to_letters,sp_str_ref)


                
                tr=tr.parse(tr_str_ref)
                sp=sp.parse(sp_str_ref)

                sp_copy= copy.deepcopy(sp)

                start_time1= time.time()

                
                tr.order_gene(sp)
                tr.label_internal()
                sp.label_internal()
                tr.map_gene(sp)
                reco.setCost(sp)
                sp.isRoot=True
                tr.isRoot=True
                reco.L_cost=  2
                reco.D_cost=  2
                sp_copy.isRoot=True
                
                
                def call_reconcILS_function():
                        return reco.iterative_reconcILS(tr,sp,sp_copy,sp,[])



                reco.gene_tree= copy.deepcopy(tr)
                species_edge_list=reco.get_edges(sp)
  
                if len(tree)<50:
                        reco.reconcILS(tr,sp,sp_copy,sp,[])
                        
                
                        li =red.sp_event(sp,[])
                else:
                        print('Using Iterative Function')
                        num_threads = 40

                        timeout_duration = 1800
                        with ThreadPoolExecutor(max_workers=num_threads) as executor:
                                future = executor.submit(call_reconcILS_function)

                        try:
                                li = future.result(timeout=timeout_duration)
                        except TimeoutError:
                                print(f"The function call timed out after {timeout_duration} seconds.")
                                li=[]
                        
                        
                        print('-----------------------------------------------------Done With reconcils-----------------------------------------------------------------------------')
                
                time_it.append( time.time()-start_time1 )

                start_time2= time.time()

                gc.collect()
                if len(li)==0:

                        continue
                else:	
                        print('##############################################################################################')
                        li.reverse()
                        dic_reconcILS={}
                        for i in li:
                                if type(i)==list:
                                        if (repr(list(sorted(i[0]))),' to ',repr(list(sorted(i[1]))))  in dic_reconcILS.keys():
                                                if i[3]==True:
                                                        dic_reconcILS[(repr(list(sorted(i[0]))),' to ',repr(list(sorted({i[1]}))))]+=[i[2]]
                                                else:
                                                        dic_reconcILS[(repr(list(sorted(i[0]))),' to ',repr(list(sorted(i[1]))))]+=[i[2]]
                                        else:   
                                                if i[3]==True:
                                                        dic_reconcILS[(repr(list(sorted(i[0]))),' to ',repr(list(sorted({i[1]}))))]=[i[2]]
                                                else:
                                                        dic_reconcILS[(repr(list(sorted(i[0]))),' to ',repr(list(sorted(i[1]))))]=[i[2]]

                        print('##############################################################################################')


                        genetree = PhyloTree(tr_str_ref)
                        sptree = PhyloTree(sp_str_ref)
                        recon_tree, events = genetree.reconcile(sptree)
                        #recon_tree.show()
                        #print(recon_tree)
                        fr={}
                        for node in recon_tree.traverse(strategy="preorder"):
                                if len(node.children) >= 0:
                                        if hasattr(node,'evoltype'):
                                                if node.evoltype in ['D','L']:
                                                        node_ =sorted(node.get_species())
                                                        if node.evoltype=='L' and len(node.children) >1:
                                                                if not hasattr(node.children[1],'evoltype') or not hasattr(node.children[0],'evoltype') or (hasattr(node.children[1],'evoltype')  and node.children[1].evoltype!='L') or  (hasattr(node.children[0],'evoltype') and node.children[0].evoltype!='L'):
                                                                                        continue
                                                                else:
                                                                        if node.up:
                                                                                node_up = sorted(node.up.get_species())
                                                                                if node.evoltype=='L' and node_up==node_:
                                                                                        continue
                                                                                elif (repr(node_up),' to ',repr(node_))  in fr.keys():
                                                                                        fr[(repr(node_up),' to ',repr(node_))]+=[node.evoltype]
                                                                                else:
                                                                                        fr[(repr(node_up),' to ',repr(node_))]=[node.evoltype]
                                                                        else:
                                                                        
                                                                                if (repr(node_),' to ',repr(node_))  in fr.keys():
                                                                                        fr[(repr(node_),' to ',repr(node_))]+=[node.evoltype]
                                                                                else:
                                                                                        fr[(repr(node_),' to ',repr(node_))]=[node.evoltype]
                                                                        if hasattr(node.children[1],'evoltype'):
                                                                                node.children[1].evoltype='S'
                                                                        if hasattr(node.children[0],'evoltype'):
                                                                                node.children[0].evoltype='S'

                                                        else:                                
                                                                if node.up:
                                                                        node_up = sorted(node.up.get_species())
                                                                        if node.evoltype=='L' and node_up==node_:
                                                                                continue
                                                                        if (repr(node_up),' to ',repr(node_))  in fr.keys():
                                                                                fr[(repr(node_up),' to ',repr(node_))]+=[node.evoltype]
                                                                        else:
                                                                                fr[(repr(node_up),' to ',repr(node_))]=[node.evoltype]
                                                                else:
                                                                
                                                                        if (repr(node_),' to ',repr(node_))  in fr.keys():
                                                                                fr[(repr(node_),' to ',repr(node_))]+=[node.evoltype]
                                                                        else:
                                                                                fr[(repr(node_),' to ',repr(node_))]=[node.evoltype]





     


                        eve_rec= {'D':0,'I':0,'L':0}
                        eve_ete= {'D':0,'I':0,'L':0}



                        for i in species_edge_list:
                               

                                
                                rec_in_dic= {'D':0,'I':0,'L':0}
                                

                                if i in dic_reconcILS:
                                        dic_in= dict(Counter(dic_reconcILS[i]))
                                        for j in rec_in_dic:
                                                if j in dic_in:
                                                        rec_in_dic[j]=dic_in[j]
                                                        eve_rec[j]+=dic_in[j]
                                        dic_reconcILS[i]=rec_in_dic
                                else:
                                        dic_reconcILS[i]={'D':0,'I':0,'L':0}
                                
  

                           
                                
                                
                                ete_in_dic= {'D':0,'I':0,'L':0}

                                if i in fr:
                                        dic_in= dict(Counter(fr[i]))
                                        for j in ete_in_dic:
                                                if j in dic_in:
                                                        ete_in_dic[j]=dic_in[j]
                                                        eve_ete[j]+=dic_in[j]
                                        fr[i]=ete_in_dic
                                else:
                                        fr[i]={'D':0,'I':0,'L':0}

                                
                                
                        sp_small=red.parse(sp_str_ref)
                        sp_small.label_internal()

                        #sp=sp.parse(sp_str_ref)
                        #sp.label_internal()

                        reco.edge_to_event(sp_large,dic_reconcILS,0) 
                        reco.edge_to_event(sp_large,fr,1)


                        reco.edge_to_event(sp_small,dic_reconcILS,0) 
                        reco.edge_to_event(sp_small,fr,1)

  
                        dic['Gene_tree']+=[red.to_newick(reco.gene_tree)]
                        dic['Species_Tree']+=[red.to_newick_sp(sp)]
                        dic['Replicate']+=[ili]
                        dic['Process']+=['reconcILS']  
                        dic['Gene_tree']+=[recon_tree.write(format=9)]
                        dic['Species_Tree']+=[red.to_newick_sp(sp)]
                        dic['Replicate']+=[ili]
                        dic['Process']+=['ete3']
                        dic['Duplication']+=[eve_rec['D'],eve_ete['D']]
                        dic['NNI']+=[eve_rec['I'],eve_ete['I']]
                        dic['Loss']+=[eve_rec['L'],eve_ete['L']]

                        print('-----------------------------------------------------Done With ETE3---------------------------------------------------------------------------')
                
                        df = pd.DataFrame(dic)
                        df.to_csv('no_dup.csv', index=False)
                        dic1= {'time':time_it}
                        df_time= pd.DataFrame(dic1)
                        df_time.to_csv('timing_no_dup.csv', index=False)
                        oer= open('labeled_L_no_dup.tre','w+')
                        oer.write(red.to_newick_sp(sp_large))
                        oer.close()
                        #oer= open('labeled_L_'+str(ili)+'.tre','w+')
                        #oer.write(red.to_newick_sp(sp_small))
                        #oer.close()
                        print('-----------------------------------------------------next_replicate-----------------------------------------------------------------------------')
                        
                        time_it2.append(time.time()-start_time2)
                        dic1= {'time':time_it2}
                        df_time= pd.DataFrame(dic1)
                        df_time.to_csv('timing_no_dup_write_1.csv', index=False)

    