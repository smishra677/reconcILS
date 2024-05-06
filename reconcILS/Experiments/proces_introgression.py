import pandas as pd
from collections import Counter
import ast
import sys
sys.path.append("../")
from reconcILS import *

df = pd.read_csv('discordance.csv', converters={'events': eval}) 
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
    'Carlitosyrichta': 'A',
    'Cebuscapucinus': 'B',
    'Cercocebusatys': 'C',
    'Macacafascicularis': 'D',
    'Macacanemestrina': 'E',
    'Theropithecusgelada': 'F',
    'Papioanubis': 'G',
    'Macacamulatta': 'H',
    'Mandrillusleucophaeus': 'I',
    'Chlorocebussabaeus': 'J',
    'Colobusangolensis': 'K',
    'Piliocolobustephrosceles': 'L',
    'Rhinopithecusbieti': 'M',
    'Rhinopithecusroxellana': 'N',
    'Gorillagorilla': 'O',
    'Homosapiens': 'P',
    'Panpaniscus': 'Q',
    'Pantroglodytes': 'R',
    'Pongoabelii': 'S',
    'Nomascusleucogenys': 'T',
    'Saimiriboliviensis': 'U',
    'Aotusnancymaae': 'V',
    'Callithrixjacchus': 'W',
}

sp_string=str(open('./sp_tree_pruned.tre').read()) 




sp_string = sp_string.replace('e-', '0')



reco =reconcils()
red= readWrite.readWrite()
sp_large=red.parse_bio(sp_string)

sp_str_ref=sp_large.to_newick()

sp_str_ref = generate_species_codes(species_to_letters,sp_str_ref)


sp_large=sp_large.parse(sp_str_ref)
sp_large.label_internal()
species_edge_list=reco.get_edges(sp_large)

dic_reconcILS={}
def write_events_sp(tree):
        if tree.isLeaf:
            ev=''
            for i in tree.event_list:
                print(i)
                ev+=str(i[0])[1:-1]
                ev+='   '
            return ev+tree.taxa   
        else:
            ev=''
            for i in tree.event_list:
                print(i)
                ev+=str(i[0])[1:-1]
                ev+='   '
                
            return ev

def traverse_sp(newick,tree):
        if tree.leftChild and not tree.rightChild:
            newick = f"(,{traverse_sp(newick,tree.leftChild)}){write_events_sp(tree)}"
        elif not tree.leftChild and tree.rightChild:
            newick = f"({traverse_sp(newick,tree.rightChild)},){write_events_sp(tree)}"
        elif tree.leftChild and tree.rightChild:
            newick = f"({traverse_sp(newick,tree.rightChild)},{traverse_sp(newick,tree.leftChild)}){write_events_sp(tree)}"
        elif not tree.leftChild and not tree.rightChild :
            newick = f"{write_events_sp(tree)}"
        else:
            pass
        return newick





# Got it From StackOverflow:
#https://stackoverflow.com/questions/61117131/how-to-convert-a-binary-tree-to-a-newick-tree-using-python
# https://stackoverflow.com/questions/61117131/how-to-convert-a-binary-tree-to-a-newick-tree-using-python
def to_newick_sp(tree):
        newick = ""
        newick = traverse_sp(newick,tree)
        newick = f"{newick};"
        return newick

for i in df['branch']:
    mask = df['branch'] == i
    value= df.loc[mask, 'events']
    dic=dict(value.values[0][0])
    deno =sum(dic.values())


    li_key = list(dic.keys())
    num = dic[li_key[0]]-dic[li_key[1]]

    delta=0
    if deno>0:
        delta = num/deno
    


    if delta<0:

        dic_reconcILS[i]={'d':delta}
    else:
        dic_reconcILS[i]={'d':delta}


new_dic_reconcILS={}
for ij in species_edge_list:
    li1= str(ij)[1:-1]
    dic_key =list(dic_reconcILS.keys())
    dic_key = [ijk[1:-1] for ijk in dic_key]

    if li1 not in dic_key:
        new_dic_reconcILS[ij]={'d':0}
    else:
        new_dic_reconcILS[ij]=dic_reconcILS['['+li1+']']




reco.edge_to_event(sp_large,new_dic_reconcILS,0) 
print(to_newick_sp(sp_large))