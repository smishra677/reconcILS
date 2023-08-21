
import Tree
import re
import os
import copy

import argparse



a= Tree.Tree()

def read_trees(i,folder):
     gene_tre= open('./'+folder+'/rep_'+str(i)+'.tre')
     tr =gene_tre.read().strip().split('\n')
     gene_tre.close()
     return str(tr[0])

def write_tree(i,folder,tree):
    with open('./'+folder+'/rep_1_'+str(i)+'.tre', 'w') as f:
        f.write(tree)
    

def read_log(flag,i,dic,folder):
    o= open('./'+folder+'/rep_'+str(i)+'.log').read()
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
        ##print(dir(root))

            print(li)
            
            li.append(node.evoltype)
        ##print(root.isLeaf),
    return li


def sp_event(root,li):
 
    if root:
        sp_event(root.leftChild,li)
        sp_event(root.rightChild,li)
        #print(root.taxa),
        #print('event',root.evolve),
        #print('cost',root.cost)
        if (root.evolve=='NNI'):
            if root.cost==0:
                 li.append(root.evolve)
            else:
                li+= ['NNI' for i in range(root.cost)]
        else:
            li.append(root.evolve)
        ##print(root.isLeaf),
    return li



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


def print_label(root):
    if root:
        if root.isLeaf:
            if root.tag:
                print(root.tag.taxa)
            else:
                print(root.tag)
        else:
            if root.event:
                print(root.event.taxa)
            else:
                print(root.event)
        print(root.taxa)

        print('$$')
        print_label(root.leftChild)

        print_label(root.rightChild)

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

def clearcost_1(sp,left,right):
    left.refTo=list(set(sp.refTo).difference(set(right.refTo)))
    right.refTo=list(set(sp.refTo).difference(set(left.refTo)))

def find_parent_child(root,child):

    if len(root.refTo)>1:
            #root.refTo.reverse()
            print('ref',root.refTo)
            for tre in root.refTo:
                for tree1 in root.refTo:
                    if tree1.parent:
                        print(tre.children,tree1,tree1.taxa,tree1.parent,tree1.parent.taxa,tree1.parent.children)
                    if (tree1 in tre.children):
                        if tree1==tre.leftChild:
                            print('match_left')
                            child.append([tre,tree1,'Left'])
                        else:
                            child.append([tre,tree1,'Right'])

    return child
        
def clear_ref(root):
    if root:
        clear_ref(root.leftChild)

        clear_ref(root.rightChild)
        root.refTo=[]
        root.tag = None
            
def parent_child(root,child):
    if root:
        if root.isLeaf:
            return []
        else:
            child= find_parent_child(root,child)
    return child


def clearid(sp,ori):
    if sp:
        if ori=='Left':
            sp.id =sp.id*5


        else:
            sp.id =sp.id*7

        clearid(sp.leftChild,ori)
            
        clearid(sp.rightChild,ori) 
            


def paralogy_NNI(sp,lis_paralogy,lis_NNI):
    if sp:
        list_paralogy,lis_NNI=paralogy_NNI(sp.leftChild,lis_paralogy,lis_NNI)
        list_paralogy,lis_NNI=paralogy_NNI(sp.rightChild,lis_paralogy,lis_NNI)
        if len(sp.paralogy)>0:
            lis_paralogy+=sp.paralogy
        if len(sp.NNI_)>0:
            lis_NNI+=sp.NNI_[0]
        return list_paralogy,lis_NNI

    return lis_paralogy,lis_NNI

            



def clearcost(sp):
    if sp:
        if sp.isLeaf:
            sp.refTo =sp.refTo[1:]

        clearcost(sp.leftChild)
        
        clearcost(sp.rightChild)  


def copy_event(re,sp):
    if sp:
        sp.NNI_.append(re.NNI_)
        copy_event(re.leftChild,sp.leftChild)
        copy_event(re.rightChild,sp.rightChild)


def sp_id(sp,dic):
    if sp:
        dic=sp_id(sp.leftChild,dic)
        dic=sp_id(sp.rightChild,dic)

        #dic[sp.id]=set(sp.taxa)
    return dic



def setCost(sp):
    if sp:
        sp.inital_ref= len(sp.refTo)
    

        setCost(sp.leftChild)
        
        setCost(sp.rightChild)  


def swap(tree):
    new_tree= tree.leftChild
    tree.leftChild= tree.righChild
    tree.rightChild= new_tree
    return tree

def ILS(gene_tree,tr,sp_copy,cost):
    list_tree=[]
    child=[]

    child= parent_child(tr,child)
    print(to_newick(gene_tree))
    print(child)
    #print(gene_tree.parent)

    if len(child)==0 or cost==0:
        return gene_tree,cost
    else:
        #child =[child[0]]
        new_topo=copy.deepcopy(gene_tree)
        geneTree =copy.deepcopy(new_topo)
        geneTree.reset()
        
        
        for ch1 in child:
                if ch1 in tr.visited:
                    print('visited')
                    continue 
                ch = copy.deepcopy(ch1[0])
                ch.reset()
                
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
                    

                    
                    #new_cost =cop.optimize_cost(i[0],i[1])
                    i[0].order_gene(cop)
                    i[0].label_internal()
                    cop.label_internal()
                    i[0].map_gene(cop)
                    new_cost = len(cop.refTo)
                    print(best_cost,new_cost)

                
                    
                    if best_cost>new_cost and cost>0:
                        best_cost=new_cost


                        if tr.isRoot:
                             
                            new_topo=copy.deepcopy(i[1])
                        else:
                             new_topo=copy.deepcopy(i[1])
                        
                        if ch1[2]=='Left':
                            if new_topo.leftChild.taxa ==None:
                                tr.NNI_+=[tr.id,tr.id]
                            else:
                                if len(new_topo.leftChild.taxa)==1:
                                    if tr.leftChild.rightChild:
                                        tr.NNI_+=[(tr.id,tr.leftChild.rightChild.id)]
                                    else:
                                        tr.NNI_+=[(tr.id,tr.rightChild.id)]
                                else:
                                    tr.NNI_+=[(tr.id,tr.leftChild.id)]
                        else:
                            if new_topo.rightChild.taxa==None:
                                tr.NNI_+=[tr.id,tr.id]
                            else:
                                if  len(new_topo.rightChild.taxa)==1:
                                    if tr.rightChild.leftChild:
                                        #print('Taxa',tr.rightChild.rightChild.taxa)
                                        tr.NNI_+=[(tr.id,tr.rightChild.leftChild.id)]
                                    else:
                                        tr.NNI_+=[(tr.id,tr.leftChild.id)]
                                else:
                                    tr.NNI_+=[(tr.id,tr.rightChild.id)]
                            
                    #print('new_topo',to_newick(new_topo))
                        imporvement=True

                cost=cost-1
                tr.visited.append([ch1[0],ch1[1]])
                
                if cost==0 or imporvement==False:
                        return new_topo,cost
                else:
                        new_sp = copy.deepcopy(sp_copy)
                        new_sp.reset()
                        new_topo.reset()
                        new_topo.order_gene(new_sp)
                        
                        new_topo.label_internal()
                    
                        new_sp.label_internal()
                        
                        new_topo.map_gene(new_sp)

                        #print(to_newick(new_sp))

                        return ILS(new_topo,new_sp,sp_copy,cost)
                
                
                '''



                print(to_newick(new_topo))

    
            
             
               

             
                '''
    return new_topo,cost
        
def reconcILS(tr,sp,sp_copy,sp_):
    if sp: 
        #print(sp.taxa)
        #print(sp.children)

        sp_copy = copy.deepcopy(sp)

        ##print(sp.evolve)


        tr_copy_1 = copy.deepcopy(tr)
        sp_1 =copy.deepcopy(sp) 

        tr_copy_2 = copy.deepcopy(tr)
        sp_2 =copy.deepcopy(sp) 


        Initial_multiple_mapping=len(sp.refTo)
        
        print('cost',Initial_multiple_mapping)

        if sp.evolve!=None:
            sp.leftChild= reconcILS(tr,sp.leftChild,sp_copy,sp_)
            sp.rightChild=reconcILS(tr,sp.rightChild,sp_copy,sp_)
            return sp        

        if tr==None:
            if sp.isLeaf:
                if sp.inital_ref==0:
                    
                    sp.evolve='Loss'
                    sp = parse(to_newick(sp))
                    sp.evolve='Loss'
                elif Initial_multiple_mapping>=2:
                    sp.evolve='Duplication'
                    sp = parse(to_newick(sp))
                
                return sp 
            else:
                if Initial_multiple_mapping>=2:
                    sp.evolve= 'Duplication'
                else:
                    sp.evolve='Speciation'
                
                #sp.leftChild= reconcILS(None,sp.leftChild,sp_copy,sp_)
                #sp.rightChild =reconcILS(None,sp.rightChild,sp_copy,sp_)
                return sp

        print(to_newick(tr))
        print(to_newick(sp))  
        if Initial_multiple_mapping in [0,1]:
            
        
            if sp.isLeaf:
                if sp.inital_ref==0 and sp.parent.evolve!='Loss':
                    sp.evolve='Loss'

                    sp.cost=0
                    sp = parse(to_newick(sp))
                    sp.evolve='Loss'
                    return sp
                else:
                    sp.evolve='Speciation'
                    sp = parse(to_newick(sp))
                    sp.evolve='Speciation'
                    return sp
            elif (sp.leftChild.isLeaf and sp.rightChild.isLeaf) :

                    sp.evolve= 'Speciation'
                    print('species',to_newick(sp))

                    if len(set(sp.leftChild.taxa).intersection(set(tr.taxa)))==0:
                        sp.leftChild.evolve='Loss'
                        clearid(sp,'right')
                        return reconcILS(tr,sp,sp_copy,sp_)
                    
                    if len(set(sp.rightChild.taxa).intersection(set(tr.taxa)))==0:
                        sp.rightChild.evolve='Loss'
                        clearid(sp,'Left')
                        return reconcILS(tr,sp,sp_copy,sp_)
                    

                    if len(set(sp.leftChild.taxa).intersection(set(tr.leftChild.taxa)))>=len(set(sp.rightChild.taxa).intersection(set(tr.leftChild.taxa))):


                        sp.leftChild= reconcILS(tr.leftChild,sp.leftChild,sp_copy,sp_)
                        sp.rightChild =reconcILS(tr.rightChild,sp.rightChild,sp_copy,sp_)
                    else:

                        sp.leftChild= reconcILS(tr.rightChild,sp.leftChild,sp_copy,sp_)
                        sp.rightChild =reconcILS(tr.leftChild,sp.rightChild,sp_copy,sp_)                        
                    
                    sp.cost=0
                    #print('1283',to_newick(sp))
                    
                    return sp
            else:
                    sp.evolve= 'Speciation'
      
                    if len(set(sp.leftChild.taxa).intersection(set(tr.taxa)))==0:
                        print('left')
                        sp.leftChild.evolve='Loss'
                        clearid(sp,'right')
                        return reconcILS(tr,sp,sp_copy,sp_)
                    
                    if len(set(sp.rightChild.taxa).intersection(set(tr.taxa)))==0:
                        print('right')
                        sp.rightChild.evolve='Loss'
                        clearid(sp,'Left')
                        return reconcILS(tr,sp,sp_copy,sp_)
                    

                    if len(set(sp.leftChild.taxa).intersection(set(tr.leftChild.taxa)))>=len(set(sp.rightChild.taxa).intersection(set(tr.leftChild.taxa))):

                        
                        sp.leftChild= reconcILS(tr.leftChild,sp.leftChild,sp_copy,sp_)
                        sp.rightChild =reconcILS(tr.rightChild,sp.rightChild,sp_copy,sp_)
                    else:

                        sp.leftChild= reconcILS(tr.rightChild,sp.leftChild,sp_copy,sp_)
                        sp.rightChild =reconcILS(tr.leftChild,sp.rightChild,sp_copy,sp_)                        
                    
                    sp.cost=0
                    #print('1283',to_newick(sp))
                    
                    return sp
           



        else:
                

                print(5)

                new_topo,cost =(ILS(tr_copy_1,sp_1,sp_1,Initial_multiple_mapping))
                #new_topo =parse(to_newick(new_topo))


                recon_1 = copy.deepcopy(sp_1)
                recon_1.reset()
                recon_1.label_internal()


                new_gene_tree =copy.deepcopy(tr_copy_2)
                new_gene_tree.reset()


                recon=Tree.Tree()
                recon= copy.deepcopy(sp_2)
                

                recon.tag_species(new_gene_tree,tr_copy_2)

                recon_left = copy.deepcopy(sp)
                #recon_left.id =recon_left.id *10
                recon_right = copy.deepcopy(sp)
                #recon_right.id =recon_right.id *10
                new_sp = copy.deepcopy(sp)

                new_sp.leftChild=recon_left
                new_sp.rightChild=recon_right

                new_sp.reset()
                recon_left.label_internal()

                recon_right.label_internal()

                
                if recon.split_list!=None:
                    recon.split_list[0].map_gene(recon_left)
                
                    recon.split_list[0].map_recon(recon_left)
                    
                    recon.split_list[1].map_gene(recon_right)
                
                    recon.split_list[1].map_recon(recon_right)
                else:
                    tr_copy_2.map_recon(recon)


                recon.clean_up()
                recon.total_cost_()



                tr_copy_2.leftChild.total_cost_()
                tr_copy_2.leftChild.clean_up()


                tr_copy_2.rightChild.total_cost_()
                tr_copy_2.rightChild.clean_up()

                

                recon_1_cost =tr_copy_2.rightChild.cost+tr_copy_2.leftChild.cost+1+len(recon_right.refTo)+len(recon_left.refTo)
                recon_left.reset()
                recon_right.reset()
                tr_copy_2.rightChild.reset()
                tr_copy_2.leftChild.reset()
               
                new_topo.label_internal()
                new_topo.map_gene(recon_1)
                
                new_topo.map_recon(recon_1)
 
                recon_1.clean_up()
                recon_1.total_cost_()

                recon_1.cost= recon_1.cost+(Initial_multiple_mapping- cost) +len(recon_1.refTo)

                
                print(recon_1_cost)

                print(recon_1.cost)
 
 
                if  recon_1.cost<recon_1_cost and  sp.isLeaf==None:
                        #print('NNI')
                        #print(to_newick(new_topo))
                        #print(to_newick(tr))
                        sp.refTo=[]
                        
                        new_topo.reset()
                        clear_ref(sp)


                        new_topo.order_gene(sp)

                                
                        new_topo.label_internal()
                    

                                        
                        new_topo.map_gene(sp)
                        print(sp.leftChild.refTo)
                        print(sp.rightChild.refTo)
                        #new_topo.id= new_topo.id*17
                        
                        sp.cost=Initial_multiple_mapping- cost
                        sp.evolve='NNI'
                        print('NNI',to_newick(new_topo))
                        copy_event(sp_1,sp)

                        #return reconcILS(new_topo,sp,sp_copy,sp_)
                       


                        print('NNI',to_newick(new_topo))
                        if len(set(new_topo.leftChild.taxa).intersection(set(sp.leftChild.taxa)))>len(set(new_topo.rightChild.taxa).intersection(set(sp.leftChild.taxa))):

                            
                            #sp.leftChild.optimize_cost(sp.leftChild,new_topo.leftChild)
                            #sp.rightChild.optimize_cost(sp.rightChild,new_topo.rightChild)

                            sp.leftChild = reconcILS(new_topo.leftChild,sp.leftChild,sp_copy,sp_) 
                            sp.rightChild = reconcILS(new_topo.rightChild,sp.rightChild,sp_copy,sp_)
                        else:
                            #sp.leftChild.optimize_cost(sp.leftChild,new_topo.rightChild)
                            #sp.rightChild.optimize_cost(sp.rightChild,new_topo.leftChild)

                            sp.leftChild = reconcILS(new_topo.rightChild,sp.leftChild,sp_copy,sp_) 
                            sp.rightChild = reconcILS(new_topo.leftChild,sp.rightChild,sp_copy,sp_)
      
                        #print('NNI',to_newick(sp))
                    
                        return sp

                if  recon_1.cost>=recon_1_cost or  sp.isLeaf:

                        print('Duplication')
                        
                        recon_left = copy.deepcopy(sp)
                        clearid(recon_left,'Left')
                        
                        
                        
                        recon_right = copy.deepcopy(sp)
                        clearid(recon_right,'right')
                        
                        
                        recon_left.reset()
                        recon_right.reset()

                        recon_left.label_internal()
                        recon_right.label_internal()


                        sp.taxa=''



                        #recon_left.isLeaf=True


                        #recon_right.isLeaf=True
                        sp.evolve='Duplication'
                        #sp.refTo=sp.refTo[1:]
       
                        if 1==1:
                            print(1)

                            
    
                            


                            #tr.leftChild.reset()
                            #tr.rightChild.reset()
                            #recon_left.isLeaf=True
                            #recon_right.isLeaf=True
                            #clearcost(recon_left)


                            #clearcost(recon_right)

                            tr.reset()


                            clear_ref(sp)

                            #sp.label_internal()

                        

                            #tr.map_gene(sp)
                            #tr.map_gene(recon_right)

                            li =sp_event(sp,[])
                            print(li)

                            print(to_newick(sp))

                            print(to_newick(tr))

                            print(to_newick(recon_left))

                            print(to_newick(recon_right))
                            
  

                            #clearcost_1(sp_copy,recon_left,recon_right)
                            print(to_newick(tr.leftChild))

                            print(to_newick(tr.rightChild))
                            print(tr.isLeaf)
                            if tr.isLeaf:
                                tr.order_gene(recon_left)
                                tr.label_internal()
   
                                tr.map_gene(recon_left)
                                sp.leftChild = reconcILS(tr,recon_left,sp_copy,sp_) 
                                tr.reset()
                                tr.order_gene(recon_right)
                                tr.label_internal()
   
                                tr.map_gene(recon_right)
                                sp.rightChild = reconcILS(tr,recon_right,sp_copy,sp_)
                                sp.children+=[sp.leftChild ]
                                sp.children+=[sp.rightChild]
                                sp.leftChild.parent=sp
                                sp.rightChild.parent=sp
                                if sp.parent==None:
                                    sp.paralogy+=[(-367,sp.id)]
                                else:
                                    if sp.isLeaf:
                                        sp.paralogy+=[(sp.id, sp.id)]  
                                    else:  
                                        sp.paralogy+=[(sp.parent.id, sp.id)]  
                                    sp.paralogy+[(sp.parent.id, sp.id)]
                                sp.isLeaf=None
                            


                            #print(sp.leftChild.refTo)
                            else:
                                tr.leftChild.order_gene(recon_left)
                                tr.rightChild.order_gene(recon_right)
                                tr.label_internal()

                                tr.leftChild.map_gene(recon_left)
                                sp.leftChild = reconcILS(tr.leftChild,recon_left,sp_copy,sp_) 
                                tr.rightChild.map_gene(recon_right)
                                sp.rightChild = reconcILS(tr.rightChild,recon_right,sp_copy,sp_)
                                sp.children+=[sp.leftChild ]
                                sp.children+=[sp.rightChild]
                                sp.leftChild.parent=sp
                                sp.rightChild.parent=sp

                                if sp.parent==None:
                                    sp.paralogy+=[(-367,sp.id)]
                                else:
                                    if sp.isLeaf:
                                        sp.paralogy+=[(sp.id, sp.id)]  
                                    else:  
                                        sp.paralogy+=[(sp.parent.id, sp.id)]
                                    #sp.paralogy+=[(sp.parent.id, sp.id)]
                                sp.isLeaf=None
                            
                           
                            #sp.leftChild=recon_left
                            #sp.rightChild=recon_right
                        
                            return sp
                        

def search_per(sp,recon):
    if sp:
        if sp.taxa==recon.taxa:
            sp.li.append(recon.evolve)
        search_per(sp.leftChild,recon)
        search_per(sp.rightChild,recon)

def per(sp,recon):
    if recon:
        search_per(sp,recon)

        per(sp,recon.leftChild)
        per(sp,recon.rightChild)
        


def make_table(lis_paralogy,lis_NNI,dic,sp_1):
    #dic=sp_id(sp_1,dic)
    print('###############################')
    for i in lis_paralogy:
        if(i[0]  in dic.keys() and i[1] in dic.keys()):
            print('Duplication on edge between:',dic[i[0]],'---->',dic[i[1]])
        else:
            key_0=-1
            key_1=-1
            for k in dic.keys():
                if k%i[0]==0 and key_0==-1:
                    key_0=k
                if k%i[1]==0 and key_1==-1:
                    key_1=k
            print('Duplication on edge between:',dic[key_0],'---->',dic[key_1])
                
            

    for i in lis_NNI:
        if(i[0]  in dic.keys() and i[1] in dic.keys()):
            print('ILS on edge between:',dic[i[0]],'------->',dic[i[1]])
        else:
                key_0=-1
                key_1=-1
                for k in dic.keys():
                    if k%i[0]==0 and key_0==-1:
                        key_0=k
                    if k%i[1]==0 and key_1==-1:
                        key_1=k
                print('ILS on edge between:',dic[key_0],'------->',dic[key_1])
                    
                
        
    print('###############################')




def parse1():
    parser = argparse.ArgumentParser(description="reconcile Gene Tree with Species Tree")
    parser.add_argument('--spTree', type=str, help="Species_Tree")
    parser.add_argument('--gTree', type=str, help="gene_Tree")
    parser.add_argument('--output', type=str, help="Location and name to output")
    args= parser.parse_args()
    return(args)




def main():
    parser = parse1()
    from collections import Counter
    import pandas as pd
    import igraph as ig
    import matplotlib.pyplot as plt
    


    sp_string=parser.spTree

    gene_tree=parser.gTree

    tr=parse(gene_tree)

    tr=parse(to_newick(tr))


    sp=parse(sp_string)
    
    sp_copy= copy.deepcopy(sp)

    sp_copy_1= copy.deepcopy(sp)
    sp_copy.reset()

    tr.order_gene(sp)

    tr.label_internal()
    sp.label_internal()
    tr.map_gene(sp)


    print('-----------------')


    #
    print(to_newick(tr))

    setCost(sp)
    sp.isRoot=True
    tr.isRoot=True
    sp_copy.isRoot=True
    reconcILS(tr,sp,sp_copy,sp)
    #print('######################33')
    
    
    li =sp_event(sp,[])
           
    sp_1=parse(sp_string)
    sp_1.label_internal()
    sp.label_internal()

    per(sp_1,sp)
    sp.map_gene(sp_1)

    

    #sp.viz()
    edges, node_,taxa_, color_,LC_=sp.find_all_edges([],[],[],[],[])
    lis_paralogy, lis_NNI= paralogy_NNI(sp,[],[])

    print(lis_paralogy)
    print(lis_NNI)
    print(len(node_))
    print(node_)
    print(edges)
    print(color_)
    print(LC_)
    g = ig.Graph(edges=edges)
    print(g)
    #g.vs["evolve"]=color_


    new_dic_LC= dict(zip(node_,LC_))
   
    Keys = list(new_dic_LC.keys())
    Keys.sort()
    sorted_dict_LC = {i: new_dic_LC[i] for i in Keys}
    print(sorted_dict_LC)
    


    '''   
    if len(lis_paralogy)>0:
        for edge  in lis_paralogy:
            g.add_edge(edge[0],edge[1],label='Duplication')
    if len(lis_NNI)>0:
        for edge  in lis_NNI:
            g.add_edge(edge[0],edge[1],label='NNI',edge_color='RED')

    '''
    
    to_delete = [i for i in range(0,max(node_)) if i  not in node_]
    g.delete_vertices(to_delete)

    new_dic= dict(zip(node_,color_))
    print(new_dic)
    
    Keys = list(new_dic.keys())
    Keys.sort()
    sorted_dict = {i: new_dic[i] for i in Keys}
    
    color_=list(sorted_dict.values())


    color_dict = {"Duplication": "blue", "Loss": "pink",'NNI':"red","Speciation":"yellow"}
    g.vs["color"] = [color_dict[evolve] for evolve in color_]
   
    
    print('taxa',taxa_)
    new_dic= dict(zip(node_,taxa_))
    print(new_dic)
    
    #make_table(lis_paralogy,lis_NNI,sorted_dict_LC,sp_copy_1)
    Keys = list(new_dic.keys())
    Keys.sort()
    sorted_dict = {i: new_dic[i] for i in Keys}
    
    g.vs["label"] = list(sorted_dict.values())
    print(list(sorted_dict.values()))
    s=str(g)



    fig, ax = plt.subplots()

    layout = g.layout_reingold_tilford(root=[int(s.split('\n')[3][0].split('--')[0])])
    ig.plot(g, layout=layout, target=ax)
    
    plt.show()


    sp.reset()
    print(to_newick(sp))
    tr.order_gene(sp)
    
    print(to_newick(sp))

    #print(Counter(li))
    dic={'Process':[],'Replicate':[],'Gene_tree':[],'Species_Tree':[],'Duplication':[],'NNI':[],'Loss':[]}
    
    dic['Gene_tree']+=[to_newick(tr)]
    dic['Species_Tree']+=[sp_string]
        
    dic= Create_pd('reconcILS',0,dict(Counter(li)),dic)
    
    df = pd.DataFrame(dic)
    print(dic)
    df.to_csv(parser.output, index=False)


if __name__ == "__main__":
    main()