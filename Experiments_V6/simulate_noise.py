
from ete3 import PhyloTree
from collections import Counter
import sys
from pathlib import Path
import argparse
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from reconcILS import *
from utils_reconcILS import *
import pandas as pd 
import numpy as np 


reco =reconcils()
red= readWrite.readWrite()
np.random.seed(42)

def enumerate_branch(tr):
    returni ={}
    stack=[tr]
    curr_index=0
    while stack:
        
        current_node=stack.pop()
        if current_node:
            if current_node.isRoot ==None and current_node.isLeaf ==None:
                if current_node.parent.leftChild ==current_node:

                    returni[curr_index]=[current_node,current_node.parent,'Left']
                else:
                    returni[curr_index]=[current_node,current_node.parent,'Right']
                curr_index=curr_index+1
            
            stack.append(current_node.leftChild)
            stack.append(current_node.rightChild)
    return returni

def parse1():
    parser = argparse.ArgumentParser(description="Add noise ")
    parser.add_argument('--file', type=str, help="Path to gene tree folder")
    parser.add_argument('--error', type=str, help="Error")
    args= parser.parse_args()
    return(args)

parser = parse1()
file___ =parser.file
error = parser.error

dupcoal_location=Path(file___)



poi = np.random.poisson(lam=float(error), size=1000)
op_tree = (dupcoal_location / 'trees.tre').read_text().splitlines()
errored_gene_trees= []
with (dupcoal_location / ("error_trees_" + str(error) + ".tre")).open("w") as file:
    for idx_,ij in enumerate(poi):
        #rando_branch = np.random.randint(ij,size=20)
        tree= str(op_tree[idx_])
        tree  = tree.replace('e-', '0')
        tr = red.parse(tree)
        tr.isRoot=True
        tr.label_internal()
        enum= enumerate_branch(tr)
        print(len(enum))
        if ij==0 or len(enum)==0 or ij>len(enum):
            errored_gene_trees.append(tr.to_newick())
            file.write(f"{tr.to_newick()}\n")
            continue
        else:
            tosses = np.random.binomial(n=1, p=ij/len(enum), size=len(enum))
            for idx,toss in enumerate(tosses):
                if toss==1:
                    get_branch =enum[idx]
                    val =get_branch[1].NNI(tr,flag=get_branch[2])
                    
                    toss_2 = np.random.binomial(n=1, p=1/2, size=1)
                    if toss_2==1:
                        tr = red.parse(val[1][1].to_newick())
                    else:
                        tr = red.parse(val[0][1].to_newick())
                    tr.isRoot=True
                    tr.label_internal()
                    enum= enumerate_branch(tr)
            

            errored_gene_trees.append(tr.to_newick())
            file.write(f"{tr.to_newick()}\n")




