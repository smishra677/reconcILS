
import Tree
import re
import os


from ete3 import TreeNode
from  ete3 import phylo
from ete3 import PhyloTree

from ete3 import NCBITaxa

a= Tree.Tree()
import dendropy


def read_trees(i):
     gene_tre= open('./output_gene/rep_'+str(i)+'.tre')
     tr =gene_tre.read()
     gene_tre.close()
     return tr

     


def tag(root):
 
    if root:
        tag(root.leftChild)
        print(root.taxa),
        print(root.event),
        if root.event !=None:
            print(root.event.refTo)
            print(root.event.taxa)
        #print(root.isLeaf),
        tag(root.rightChild)



def evolve(root):
 
    if root:
        evolve(root.leftChild)
        print(root.taxa),
        print(root.evolve),
        print(root.tag),
        #print(root.isLeaf),
        evolve(root.rightChild)


def sp_tag(root):
 
    if root:
        sp_tag(root.leftChild)
        print(root.taxa),
        print(root.refTo),
        #print(root.isLeaf),
        sp_tag(root.rightChild)

import copy

def sp_event(root,li):
 
    if root:
        sp_event(root.leftChild,li)
        sp_event(root.rightChild,li)
        print(root.taxa),
        print('event',root.evolve),
        print('cost',root.cost)
        li.append(root.evolve)
        #print(root.isLeaf),
    return li


def printInorder(root):
 
    if root:
        printInorder(root.leftChild)
        print(root.taxa),
        print(root.isLeaf),
        printInorder(root.rightChild)

def printorder(root):
 
    if root:
            print(root.tag),
            print(root.taxa),
            print(root.id),
            print(root.parent),
            print(root.refTo),
            print('###################'),
            printorder(root.leftChild)
            printorder(root.rightChild)
 

def to_newick(tree):
    newick = ""
    newick = traverse(tree, newick)
    newick = f"{newick};"
    return newick

def traverse(tree, newick):
    if tree.leftChild and not tree.rightChild:
        newick = f"(,{traverse(tree.leftChild, newick)}){tree.taxa if tree.isLeaf else ''}"
    elif not tree.leftChild and tree.rightChild:
        newick = f"({traverse(tree.rightChild, newick)},){tree.taxa if tree.isLeaf else ''}"
    elif tree.leftChild and tree.rightChild:
        newick = f"({traverse(tree.rightChild, newick)},{traverse(tree.leftChild, newick)}){tree.taxa if tree.isLeaf else ''}"
    elif not tree.leftChild and not tree.rightChild :
        newick = f"{tree.taxa if tree.isLeaf else ''}"
    else:
        pass
    return newick


def parse(newick):


    tokens = re.finditer(r"([^:;,()\s]*)(?:\s*:\s*([\d.]+)\s*)?([,);])|(\S)", newick+";")

    def recurse(tre,nextid = 0, parentid = -1): # one node
        thisid = nextid
        #tre.id= nextid
        children = []

        name, length, delim, ch = next(tokens).groups(0)
        tre.taxa= name
        if name!=0:
            tre.taxa= name.split('_')[0]
            tre.isLeaf= True
        if ch == "(":
            while ch in "(,":
                new_tre= Tree.Tree()
                node, ch, nextid,tre1 = recurse(new_tre,nextid+1, thisid)
                children.append(node)
                tre.children.append(tre1)



            if len(tre.children)==1:
                 tre.children =[]
                 tre =tre1
            else:
                tre.leftChild= tre.children[0]
                tre.children[0].parent =tre
                tre.children[1].parent =tre
                tre.rightChild=tre.children[1]
            



                
            name, length, delim, ch = next(tokens).groups(0)
        return {"id": thisid, "name": name, "length": length if length else None, 
                "parentid": parentid, "children": children}, delim, nextid,tre

    tre= Tree.Tree()
    val =recurse(tre)

    return val[-1]


def find_parent_child(root,child):

    if len(root.refTo)>1:
            for tre in root.refTo:
                for tree1 in root.refTo:
                    if (tree1 in tre.children):
                        if tree1==tre.leftChild:
                            child.append([tre,tree1,'Left'])
                        else:
                            child.append([tre,tree1,'Right'])

    return child
        

            
def parent_child(root,child):
    if root:
        if root.isLeaf:
            return []
        else:
            child= find_parent_child(root,child)
        #child=find_parent_child(root.leftChild,child)
        #child= find_parent_child(root.rightChild,child)
    return child


'''
new_gene_tree =copy.deepcopy(tr)
new_gene_tree.reset()
recon= copy.deepcopy(sp)
sp_tag(recon)
print('#########3')
recon.tag_species(new_gene_tree,tr)

print(to_newick(recon))
if recon.split_list!=None:
    for i in recon.split_list:
        tr.map_recon(i)
else:
    tr.map_recon(recon)

print(to_newick(recon))

recon.clean_up()
recon.total_cost_()




print(sp.find_cost(tr,0))

print('Gene_tree',to_newick(tr))
'''
def clearcost(sp):
    if sp:
        sp.refTo =sp.refTo[1:]

        clearcost(sp.leftChild)
        
        clearcost(sp.rightChild)  



def setCost(sp):
    if sp:
        sp.inital_ref= len(sp.refTo)
    

        setCost(sp.leftChild)
        
        setCost(sp.rightChild)  




def ILS(gene_tree,tr,sp_copy,cost):
    list_tree=[]
    child=[]
    
    
    
    child= parent_child(tr,child)

    print(gene_tree.parent)
    if len(child)==0 or cost==0:
        return gene_tree,cost
    else:
        #child =[child[0]]
        for ch1 in child:
                ch = copy.deepcopy(ch1[0])
                ch.reset()

                geneTree =copy.deepcopy(gene_tree)
                geneTree.reset()
                list_tree= ch1[0].NNI(geneTree,ch1[2])
                best_cost=cost
                imporvement=False
                new_topo=copy.deepcopy(geneTree)
                for i in list_tree:
                    i[1].reset()
                    i[0].reset()
                    cop= copy.deepcopy(sp_copy)
                    cop.reset()
                        
                        
                    print('top_0',to_newick(i[0]))
                    print('topo_1',to_newick(i[1]))
                    print('Species',to_newick(cop))
                    print('##############################')
                    
                    new_cost =cop.optimize_cost(i[0],i[1])
                    print(best_cost)
                    print(cost)
                    print(new_cost)
                    #print(cost)
                    #print(new_cost,best_cost)
                    if best_cost>new_cost and cost>0:
                        best_cost=new_cost
                        if tr.parent==None:
                             
                            new_topo=copy.deepcopy(i[1])
                        else:
                             new_topo=copy.deepcopy(i[0])
                        #print('new_topo',to_newick(new_topo))
                        imporvement=True

                if imporvement:
                        print('new_topo1',to_newick(new_topo))

                        cost=cost-1
                        if cost==0:
                            return new_topo,cost
                        else:
                            new_sp = copy.deepcopy(sp_copy)
                            new_sp.reset()
                            new_topo.reset()
                            new_topo.printorder_gene(new_sp)
                        
                            new_topo.label_internal()
                    
                            new_sp.label_internal()
                        
                            new_topo.map_gene(new_sp)

                            return ILS(new_topo,new_sp,sp_copy,cost)

                        #print(cop.find_cost(tr))
                else:
                    return new_topo,cost
        
def driver(tr,sp,sp_copy,sp_):
    if sp: 
        print(sp.taxa)
        print(sp.parent)
        print(to_newick(tr))
        print(sp.evolve)


        tr_copy_1 = copy.deepcopy(tr)
        sp_1 =copy.deepcopy(sp) 

        tr_copy_2 = copy.deepcopy(tr)
        sp_2 =copy.deepcopy(sp) 
        
        if tr==None:
            print('x')
            sp.evolve= 'Speciation'
            return sp


        initial_cost=len(sp.refTo)
        print('cost',initial_cost)

        if initial_cost==2:
            
            if sp.isLeaf:
                print(12455)
                if sp.parent.evolve!='Duplication':
                    print(12455)
                    sp.evolve='Duplication'
                else:
                     sp.evolve='Speciation'
                return sp
            else:
                 
                new_topo,cost =(ILS(tr_copy_1,sp_1,sp_1,initial_cost))
                recon_1 = copy.deepcopy(sp_1)

                

                new_gene_tree =copy.deepcopy(tr_copy_2)
                new_gene_tree.reset()


                recon=Tree.Tree()
                recon= copy.deepcopy(sp_2)
                

                recon.tag_species(new_gene_tree,tr_copy_2)


    
                if recon.split_list!=None:
                    for i in recon.split_list:
                        tr_copy_2.map_recon(i)
                else:
                    tr_copy_2.map_recon(recon)


                #recon.clean_up()
                recon.total_cost_()
               
    
               

                

               
                new_topo.label_internal()
                new_topo.map_recon(recon_1)
                #recon_1.clean_up()
                recon_1.total_cost_()
                recon_1.cost= recon_1.cost+initial_cost-(initial_cost-cost)

                recon_left = copy.deepcopy(sp)
                recon_right = copy.deepcopy(sp)
                new_sp = copy.deepcopy(sp)
                new_sp.leftChild=recon_left
                new_sp.rightChild=recon_right
                new_sp.reset()
                tr_copy_2.reset()

                recon_left.label_internal()
                tr_copy_2.label_internal()
                tr_copy_2.map_recon(recon_left)
                #recon_left.clean_up()
                recon_left.total_cost_()
                tr_copy_2.reset()

                tr_copy_2.label_internal()
                recon_right.label_internal()
                tr_copy_2.map_recon(recon_right)
                recon_right.total_cost_()
                recon_1_cost =recon_left.cost+1+recon_right.cost
                recon_left.reset()
                recon_right.reset()
                print(recon_1_cost)
                print(recon_1.cost)
                print(to_newick(recon_right))
                print(to_newick(recon_1))
                print(to_newick(new_topo))
                print(to_newick(recon_left))
                
 
                if  recon_1.cost<=recon_1_cost:
                        print('NNI')
                        print(to_newick(new_topo))
                        print(to_newick(tr))
                        sp.refTo=[]
                        sp.evolve='NNI'

    
                        sp.leftChild = driver(new_topo.leftChild,recon_1.leftChild,sp_copy,sp_) 
                        sp.rightChild = driver(new_topo.rightChild,recon_1.rightChild,sp_copy,sp_)
                        
                        sp.cost=1
                        
                        return sp

                if  recon_1.cost>recon_1_cost:
                        print('Duplication')
                        sp.refTo=[]

                        recon_left = copy.deepcopy(sp)
                        recon_right = copy.deepcopy(sp)
                        recon_left.reset()
                        recon_right.reset()
                        tr.leftChild.reset()
                        tr.rightChild.reset()

                        recon_left.parent=sp
                        recon_right.parent=sp

                        sp.evolve='Duplication'
                        recon_left.parent=sp
                        recon_right.parent=sp


                        print(to_newick(tr.leftChild))
                        recon_left.optimize_cost(recon_left,tr.leftChild)
                        recon_right.optimize_cost(recon_right,tr.rightChild)


                        recon_left.refTo=recon_left.refTo[2:]

                        recon_right.refTo=recon_right.refTo[2:]

                        clearcost(recon_left)

                        clearcost(recon_right)
                        
                        sp.leftChild = driver(tr.leftChild,recon_left,sp_copy,sp_) 
                        sp.rightChild = driver(tr.rightChild,recon_right,sp_copy,sp_)
                        recon_left.evolve='Speciation'
                        recon_right.evolve='Speciation'
                        sp.cost=1
                        return sp

        
            
           
        if initial_cost in [0,1]:
            
            if sp.isLeaf:
                if len(sp.refTo)!=0 :
                    print(1)
                    sp.evolve='Speciation'
                    sp.cost=0
                    return sp
                else:
                    print(23)
                    if sp.parent.evolve not in ['Duplication'] and sp.inital_ref == initial_cost:
                        sp.evolve='Loss'
                        sp.cost=1
                        return sp
                    else:
                         sp.evolve='Speciation'
                         return sp
            elif (sp.leftChild.isLeaf and sp.rightChild.isLeaf)  and len(sp.refTo)!=0 :

                sp.evolve= 'Speciation'
                sp.leftChild= driver(tr,sp.leftChild,sp_copy,sp_)
                sp.rightChild =driver(tr,sp.rightChild,sp_copy,sp_)
                return sp
            else:
                print(2)
                
                if len(sp.refTo)==0:
                    print(to_newick(sp))
                    sp.evolve='Speciation'
                    sp.leftChild= driver(tr,sp.leftChild,sp_copy,sp_)
                    sp.rightChild =driver(tr,sp.rightChild,sp_copy,sp_)
                    sp.cost=0
                    
                    return sp

                else:
                    print('1283')
                    sp.evolve='Speciation'
                    sp.leftChild= driver(tr,sp.leftChild,sp_copy,sp_)
                    sp.rightChild=driver(tr,sp.rightChild,sp_copy,sp_)
                    
                    sp.cost=0

                    return sp
           



        else:
            
            if sp.isLeaf:
                        print(4)
                        print(48)
                        print('Duplication')
                    
                        sp.refTo=[]

                        sp.evolve='Duplication'
    
                        recon_left = copy.deepcopy(sp)
                        recon_right = copy.deepcopy(sp)
                        recon_left.reset()
                        recon_right.reset()
                        tr.leftChild.reset()
                        tr.rightChild.reset()
                        
                        recon_left.optimize_cost(recon_left,tr.leftChild)
                        recon_right.optimize_cost(recon_right,tr.rightChild)
                        
                        sp.isLeaf=None
                        

                        recon_left.parent=sp
                        recon_right.parent=sp

                        recon_left.refTo=recon_left.refTo[2:]

                        recon_right.refTo=recon_right.refTo[2:]

                        clearcost(recon_left)

                        clearcost(recon_right)


                        sp.leftChild = driver(tr.leftChild,recon_left,sp_copy,sp_) 
                        sp.rightChild = driver(tr.rightChild,recon_right,sp_copy,sp_)

   


                        sp.cost=1
                        return sp
            elif (sp.leftChild.isLeaf and sp.rightChild.isLeaf):
                        print(44)
                        print('Duplication')
                        sp.refTo=[]
     
                        recon_left = copy.deepcopy(sp)
                        recon_right = copy.deepcopy(sp)
                        recon_left.reset()
                        recon_right.reset()
                        tr.leftChild.reset()
                        tr.rightChild.reset()
                        recon_left.optimize_cost(recon_left,tr.leftChild)
                        recon_right.optimize_cost(recon_right,tr.rightChild)

                        recon_left.parent=sp
                        recon_right.parent=sp

                        recon_left.refTo=recon_left.refTo[2:]

                        recon_right.refTo=recon_right.refTo[2:]


                        clearcost(recon_left)

                        clearcost(recon_right)

                        sp.evolve='Duplication'

                        sp.leftChild = driver(tr.leftChild,recon_left,sp_copy,sp_) 
                        sp.rightChild = driver(tr.rightChild,recon_right,sp_copy,sp_)
                        


                        
                        sp.cost=1
                        return sp

            else:
                print(5)
                new_topo,cost =(ILS(tr_copy_1,sp_1,sp_1,initial_cost))
                recon_1 = copy.deepcopy(sp_1)

                

                new_gene_tree =copy.deepcopy(tr_copy_2)
                new_gene_tree.reset()


                recon=Tree.Tree()
                recon= copy.deepcopy(sp_2)
                

                recon.tag_species(new_gene_tree,tr_copy_2)


    
                if recon.split_list!=None:
                    for i in recon.split_list:
                        tr_copy_2.map_recon(i)
                else:
                    tr_copy_2.map_recon(recon)


                #recon.clean_up()
                recon.total_cost_()
    
               

                

               
                new_topo.label_internal()
                new_topo.map_recon(recon_1)
                #recon_1.clean_up()
                recon_1.total_cost_()
                recon_1.cost= recon_1.cost+initial_cost-(initial_cost-cost)
                
                recon_left = copy.deepcopy(sp)
                recon_right = copy.deepcopy(sp)
                new_sp = copy.deepcopy(sp)
                new_sp.leftChild=recon_left
                new_sp.rightChild=recon_right
                new_sp.reset()
                tr_copy_2.reset()

                recon_left.label_internal()
                tr_copy_2.label_internal()
                tr_copy_2.map_recon(recon_left)
                #recon_left.clean_up()
                recon_left.total_cost_()
                tr_copy_2.reset()

                tr_copy_2.label_internal()
                recon_right.label_internal()
                tr_copy_2.map_recon(recon_right)
                recon_right.total_cost_()
                recon_1_cost =recon_right.cost +recon_left.cost
                recon_left.reset()
                recon_right.reset()
                print(recon_1_cost)
                print(recon_1.cost)
 
                if  recon_1.cost<=recon_1_cost:
                        print('NNI')
                        print(to_newick(new_topo))
                        print(to_newick(tr))
                        sp.refTo=[]
                        sp.evolve='NNI'
                        sp.leftChild = driver(new_topo.leftChild,recon_1.leftChild,sp_copy,sp_) 
                        sp.rightChild = driver(new_topo.rightChild,recon_1.rightChild,sp_copy,sp_)
                        
                        sp.cost=1
                        
                        return sp

                if  recon_1.cost>recon_1_cost:
                        print('Duplication')
                        sp.refTo=[]

                        recon_left = copy.deepcopy(sp)
                        recon_right = copy.deepcopy(sp)
                        recon_left.reset()
                        recon_right.reset()
                        tr.leftChild.reset()
                        tr.rightChild.reset()

                        recon_left.parent=sp
                        recon_right.parent=sp

                        sp.evolve='Duplication'

                        recon_left.optimize_cost(recon_left,tr.leftChild)
                        recon_right.optimize_cost(recon_right,tr.rightChild)

                        recon_left.refTo=recon_left.refTo[2:]

                        recon_right.refTo=recon_right.refTo[2:]

                        clearcost(recon_left)

                        clearcost(recon_right)

                        sp.leftChild = driver(tr.leftChild,recon_left,sp_copy,sp_) 
                        sp.rightChild = driver(tr.rightChild,recon_right,sp_copy,sp_)

                        recon_left.evolve='Speciation'
                        recon_right.evolve='Speciation'
                        
                        sp.cost=1
                        return sp

            sp.evolve='Speciation'
            sp.leftChild =driver(tr,sp.leftChild,sp_copy,sp_)
            sp.rightChild = driver(tr,sp.rightChild,sp_copy,sp_)
            return sp



import pandas as pd          
def read_log(flag,i,dic):
    o= open('./output_gene/rep_'+str(i)+'.log').read()
    dic['Process']+=[flag]
    dic['Replicate']+=[i]
    
    for i in o.split('\n'):
        val= (i.split(':'))
        if val[0]=='Total duplications':
            dic['Duplication']+=[int(val[1])]
        elif val[0]=='Total ILS':
            dic['NNI']+=[int(val[1])]
        elif val[0]=='Total losses':
            dic['Loss']+=[int(val[1])]
        else:
            continue
   
    return dic


def Create_pd(flag,i,o,dic):
    
    dic['Process']+=[flag]
    dic['Replicate']+=[i]
    dic['Duplication'].append(0)
    dic['NNI'].append(0)
    dic['Loss'].append(0)


    for i in o:
         if i in ['Duplication','NNI','Loss']:
              dic[i][-1]=o[i]
              
         
    

    return dic

def Create_pd_ete3(flag,i,o,dic):
    
    dic['Process']+=[flag]
    dic['Replicate']+=[i]
    dic['Duplication'].append(0)
    dic['NNI'].append(0)
    dic['Loss'].append(0)


    for i in o:
         if i in ['D','L']:
                if i =='D':
                   dic['Duplication'][-1]=o[i]
                elif i=='L':
                    dic['Loss'][-1]=o[i]
                else:
                    continue

    return dic
                     
            
              
         


def sp_event_ete3(root):
    li=[]
    for node in root.traverse(strategy="postorder"):
        if len(node.children) != 0:
        #print(dir(root))

            print(li)
            
            li.append(node.evoltype)
        #print(root.isLeaf),
    return li

dic={'Process':[],'Replicate':[],'Gene_tree':[],'Species_Tree':[],'Duplication':[],'NNI':[],'Loss':[]}
    
from collections import Counter




sp_string='(A,(B,C));'

tree='((A,B),(C,(B,C));'
tr=parse(tree)
sp=parse(sp_string)
sp_copy= copy.deepcopy(sp)
sp_copy.reset()

tr.printorder_gene(sp)

tr.label_internal()
sp.label_internal()



tr.map_gene(sp)

setCost(sp)
sp.isRoot=True
tr.isRoot=True
sp_copy.isRoot=True
driver(tr,sp,sp_copy,sp)
print('######################33')
print(to_newick(sp))

li =sp_event(sp,[])


exit()
#sp ='(((A,B),C),D);'
for i in range(100):
    print()



    tree=read_trees(i)
    tr=parse(tree)
    print(to_newick(tr))

    
    sp=parse(sp_string)

    dic['Gene_tree']+=[to_newick(tr)]
    dic['Species_Tree']+=[to_newick(sp)]
    
    dic= read_log('True Process',i,dic)


    dic['Gene_tree']+=[to_newick(tr)]
    dic['Species_Tree']+=[to_newick(sp)]
    sp_copy= copy.deepcopy(sp)
    sp_copy.reset()

    tr.printorder_gene(sp)

    tr.label_internal()
    sp.label_internal()



    tr.map_gene(sp)


    sp.isRoot=True
    tr.isRoot=True
    sp_copy.isRoot=True
    driver(tr,sp,sp_copy,sp)
    print('######################33')
    print(to_newick(sp))

    li =sp_event(sp,[])

    dic= Create_pd('Our_algorithm',i,dict(Counter(li)),dic)

    dic['Gene_tree']+=[to_newick(tr)]
    dic['Species_Tree']+=[sp_string]


    genetree = PhyloTree(to_newick(tr))
    sptree = PhyloTree(sp_string)

    recon_tree_ete, events = genetree.reconcile(sptree)

    dic= Create_pd_ete3('ETE3',i,dict(Counter(sp_event_ete3(recon_tree_ete))),dic)


df = pd.DataFrame(dic)

df.to_csv('result.csv', index=False)

'''

initial_cost=1
#sp.find_cost(tr,0)
new_topo,cost =(ILS(tr,sp,sp_copy,initial_cost))
print(to_newick(new_topo))
print(initial_cost-cost)




#tr.map_recon(recon)
print(to_newick(recon))
new_tr= copy.deepcopy(tr)
new_tr.reset()
print(new_tr.optimize_cost(recon,recon))
#recon.tag_loss()

print(to_newick(recon))

'''