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

sp_string=str(open('./sp.tre').read()) # species Tree string
gT_list = open('./ZeroCol_ASTRAL_ML_SCO_MIN27.tre')
gT_list_r= gT_list.read()
gT_list.close()
gT_list=gT_list_r.split(';')

sort_ = sorted(enumerate(gT_list), key=lambda x: len(x[1]))

index_mapping = {original_index: new_index for new_index, (original_index, _) in enumerate(sort_)}


sp_string = sp_string.replace('e-', '0')


#dic={'Process':[],'Replicate':[],'Gene_tree':[],'Species_Tree':[],'Duplication':[],'NNI':[],'Loss':[]}
dic={'Process':[],'Replicate':[],'Gene_tree':[],'Species_Tree':[],'Duplication':[],'NNI':[],'Loss':[]}
DATA_PATH = "./" # path to the gene tree folder
file_path = os.path.join(DATA_PATH)
time_it=[]


reco =reconcils()
red= readWrite.readWrite()
sp_large=red.parse_bio(sp_string)

sp_str_ref=sp_large.to_newick()

sp_str_ref = generate_species_codes(species_to_letters,sp_str_ref)


sp_large=sp_large.parse(sp_str_ref)
sp_large.label_internal()



import gc

for ili in range(len(sort_)):
                print('##########################################################################################################################')
                print(ili)
                print('--------------------------------------------------------------------------------------------------------------------------')

                reco =reconcils()
                red= readWrite.readWrite()
        
                ili,_=sort_[ili]
        
                tree= str(gT_list[ili])

                if len(tree)==0:
                        continue

                tree = tree.replace('e-', '0')

                
                
                tr= red.parse_bio(red.to_newick(red.parse_bio(tree)))

                sp=red.parse_bio(sp_string)
                #dic= red.read_log('True Process',i,dic,gene_folder)
                tr_str_ref=tr.to_newick()
                sp_str_ref=sp.to_newick()
                tr_str_ref = generate_species_codes(species_to_letters,tr_str_ref)
                sp_str_ref = generate_species_codes(species_to_letters,sp_str_ref)


                
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
                reco.L_cost=  1.5
                sp_copy.isRoot=True

                


                reco.gene_tree= copy.deepcopy(tr)
                species_edge_list=reco.get_edges(sp)
  
                if len(tree)<50:
                        reco.reconcILS(tr,sp,sp_copy,sp,[])
                        
                
                        li =red.sp_event(sp,[])
                else:
                        print('Using Iterative Function')
                        li= reco.iterative_reconcILS(tr,sp,sp_copy,sp,[])
                        print('-----------------------------------------------------Done With reconcils-----------------------------------------------------------------------------')
                
                time_it.append( time.time()-start_time1 )


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
                                                dic_reconcILS[(repr(list(sorted(i[0]))),' to ',repr(list(sorted(i[1]))))]+=[i[2]]
                                        else:
                                                dic_reconcILS[(repr(list(sorted(i[0]))),' to ',repr(list(sorted(i[1]))))]=[i[2]]
                        print('##############################################################################################')


                        genetree = PhyloTree(tr_str_ref)
                        sptree = PhyloTree(sp_str_ref)
                        recon_tree, events = genetree.reconcile(sptree)
                        fr={}
                        for node in recon_tree.traverse(strategy="preorder"):
                                if len(node.children) >= 0:
                                        if hasattr(node,'evoltype'):
                                                if node.evoltype in ['D','L']:
                                                        node_ =sorted(node.get_species())
                                                        if node.evoltype=='L' and len(node.children) >1:
                                                                if (hasattr(node.children[1],'evoltype')  and node.children[1]!='L') or  (hasattr(node.children[0],'evoltype') and node.children[0]!='L'):
                                                                        continue
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

                                
                                


                        #sp=sp.parse(sp_str_ref)
                        #sp.label_internal()

                        reco.edge_to_event(sp_large,dic_reconcILS,0) 
                        reco.edge_to_event(sp_large,fr,1)
                        print(red.to_newick_sp(sp_large))  
  
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
                
                        print(dic)
                        df = pd.DataFrame(dic)
                        df.to_csv('bio_result_ZeroCol_ASTRAL_ML_SCO_MIN27.csv', index=False)
                        dic1= {'time':time_it}
                        df_time= pd.DataFrame(dic1)
                        df_time.to_csv('timing_ZeroCol_ASTRAL_ML_SCO_MIN27.csv', index=False)
                        oer= open('labeled_L_ZeroCol_ASTRAL_ML_SCO_MIN27.tre','w+')
                        oer.write(red.to_newick_sp(sp_large))
                        oer.close()
                        print('-----------------------------------------------------next_replicate-----------------------------------------------------------------------------')
                
                


    