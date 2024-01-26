import os
import copy
from ete3 import PhyloTree  
from collections import Counter
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
gT_list = open('./ZeroCol_ASTRAL_ML_ALLPARALOGS_MIN27.tre')
gT_list_r= gT_list.read()
gT_list.close()
gT_list=gT_list_r.split(';')[:100]


sp_string = sp_string.replace('e-', '0')

dic={'Process':[],'Replicate':[],'Gene_tree':[],'Species_Tree':[],'Duplication':[],'NNI':[],'Loss':[]}

DATA_PATH = "./" # path to the gene tree folder
file_path = os.path.join(DATA_PATH)



for i in range(len(gT_list)):

                reco =reconcils()
                red= readWrite.readWrite()
                #dic= red.read_log('True Process',i,dic,gene_folder)

        
            
        #try:
                tree= str(gT_list[i])
                tree = tree.replace('e-', '0')

                
                
                tr= red.parse_bio(red.to_newick(red.parse_bio(tree)))

                sp=red.parse_bio(sp_string)
                #dic= red.read_log('True Process',i,dic,gene_folder)
                tr_str_ref=tr.to_newick()
                sp_str_ref=sp.to_newick()
                tr_str_ref = generate_species_codes(species_to_letters,tr_str_ref)
                sp_str_ref = generate_species_codes(species_to_letters,sp_str_ref)


                genetree = PhyloTree(tr_str_ref)

                sptree = PhyloTree(sp_str_ref)  

                recon_tree, events = genetree.reconcile(sptree)

                li1=red.sp_event_ete3(recon_tree)
                dic['Gene_tree']+=[tr.to_newick()]
                dic['Species_Tree']+=[sp_string]
                        
                
                red.Create_pd_ete3('ete3',0,Counter(li1),dic)

     

                sp_copy= copy.deepcopy(sp)

                

                tr.order_gene(sp)
                tr.label_internal()
                sp.label_internal()
                tr.map_gene(sp)
                reco.setCost(sp)
                sp.isRoot=True
                tr.isRoot=True
                sp_copy.isRoot=True

                


                reco.gene_tree= copy.deepcopy(tr)

                if len(tree)<50:
                        reco.reconcILS(tr,sp,sp_copy,sp,[])
                        
                
                        li =red.sp_event(sp,[])
                else:
                        print('Using Iterative Function')
                        li= reco.iterative_reconcILS(tr,sp,sp_copy,sp,[])
                
        
                #reco.reconcILS(tr,sp,sp_copy,sp,[])

                
                dic['Gene_tree']+=[red.to_newick(reco.gene_tree)]
                dic['Species_Tree']+=[sp_string]
                #li =red.sp_event(sp,[])

                dic= red.Create_pd('reconcILS',i,li,dic)

                df = pd.DataFrame(dic)
                print(df)


                df.to_csv('bio_result.csv', index=False)

        #except:
                #continue
    



