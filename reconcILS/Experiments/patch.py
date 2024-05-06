import pandas as pd
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

df = pd.read_csv('bio_result_ZeroCol_ASTRAL_ML_ALLPARALOGS_MIN27_1.csv')


filtered_rows = df[df['Process'] == 'ete3'].index
red= readWrite.readWrite()
sp_large=red.parse_bio('(AC,(B,(A,((F,(D,E)),(C,(((AB,AA),(G,Z)),((((P,Q),(R,S)),(O,((J,(I,M)),((H,N),(L,K))))),(Y,(X,(T,((V,W),U)))))))))));')

for index in filtered_rows:
                        print(index)
                        dic={}
                        reco =reconcils()
                        red= readWrite.readWrite()
                        Gene_tree = df.at[index-1,'Gene_tree']
                        Species_Tree= '(AC,(B,(A,((F,(D,E)),(C,(((AB,AA),(G,Z)),((((P,Q),(R,S)),(O,((J,(I,M)),((H,N),(L,K))))),(Y,(X,(T,((V,W),U)))))))))));'
                        tr= red.parse_bio(Gene_tree)

                        sp=red.parse_bio(Species_Tree)
                        #dic= red.read_log('True Process',i,dic,gene_folder)
                        tr_str_ref=tr.to_newick()
                        sp_str_ref=sp.to_newick()
                        print(tr_str_ref)
                        print(sp_str_ref)


                        
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
                        reco.gene_tree= copy.deepcopy(tr)
                        species_edge_list=reco.get_edges(sp)

                        genetree = PhyloTree(tr_str_ref)
                        sptree = PhyloTree(sp_str_ref)
                        recon_tree, events = genetree.reconcile(sptree)
                        #recon_tree.show()
                        #print(recon_tree)
                        fr={}
                        dic_reconcILS={}
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

                        
                        #print(dic_reconcILS)
                        reco.edge_to_event(sp_large,dic_reconcILS,0) 
                        reco.edge_to_event(sp_large,fr,1)  


                        df.at[index, 'Duplication'] = eve_ete['D'] 
                        df.at[index, 'Loss'] = eve_ete['L'] 
 
                        print('-----------------------------------------------------Done With ETE3---------------------------------------------------------------------------')
                

                        oer= open('labeled_L_ZeroCol_ASTRAL_ML_SCO_MIN27_1_ete.tre','w+')
                        oer.write(red.to_newick_sp(sp_large))
                        oer.close()
                       

                

    #df.at[index, 'C'] = new_value 
 

df.to_csv('bio_result_ZeroCol_ASTRAL_ML_ALLPARALOGS_MIN27_1_ete.csv', index=False)