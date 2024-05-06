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

sp_string=str(open('./sp.tre').read()) 
gT_list = open('./ZeroCol_ASTRAL_ML_ALLPARALOGS_MIN27.tre')
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

time_it2=[]

reco =reconcils()
red= readWrite.readWrite()
sp_large=red.parse_bio(sp_string)

sp_str_ref=sp_large.to_newick()

sp_str_ref = generate_species_codes(species_to_letters,sp_str_ref)


sp_large=sp_large.parse(sp_str_ref)
sp_large.label_internal()
red.write_introgression(sp_large)


import gc


for ili in range(8358,len(sort_)):
