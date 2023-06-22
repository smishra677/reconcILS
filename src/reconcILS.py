
import Tree
import re
import os
import copy



a= Tree.Tree()

def read_trees(i,folder):
     gene_tre= open('./'+folder+'/rep_'+str(i)+'.tre')
     tr =gene_tre.read().strip().split('\n')
     gene_tre.close()
     return str(tr[0])

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
    return child



def clearcost(sp):
    if sp:
        if sp.isLeaf:
            sp.refTo =sp.refTo[1:]

        clearcost(sp.leftChild)
        
        clearcost(sp.rightChild)  



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
    #print(gene_tree.parent)
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
                        
                        
                    ##print('top_0',to_newick(i[0]))
                    ##print('topo_1',to_newick(i[1]))
                    ##print('Species',to_newick(cop))
                    ##print('##############################')
                    
                    new_cost =cop.optimize_cost(i[0],i[1])
                    ##print(best_cost)
                    ##print(cost)
                    ##print(new_cost)
                    ##print(cost)
                    ##print(new_cost,best_cost)
                    if best_cost>new_cost and cost>0:
                        best_cost=new_cost
                        if tr.parent==None:
                             
                            new_topo=copy.deepcopy(i[1])
                        else:
                             new_topo=copy.deepcopy(i[0])
                        ##print('new_topo',to_newick(new_topo))
                        imporvement=True
                cost=cost-1
                if cost==0 or imporvement==False:
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

                if imporvement:
                        print('new_topo1',to_newick(new_topo))
                else:
                    return new_topo,cost
        
def driver(tr,sp,sp_copy,sp_):
    if sp: 
        ##print(sp.taxa)
        ##print(sp.parent)

        sp_copy = copy.deepcopy(sp)

        ##print(sp.evolve)


        tr_copy_1 = copy.deepcopy(tr)
        sp_1 =copy.deepcopy(sp) 

        tr_copy_2 = copy.deepcopy(tr)
        sp_2 =copy.deepcopy(sp) 


        Initial_multiple_mapping=len(sp.refTo)
        #print('cost',Initial_multiple_mapping)
        


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

                sp.leftChild= driver(None,sp.leftChild,sp_copy,sp_)
                sp.rightChild =driver(None,sp.rightChild,sp_copy,sp_)
                return sp

          
        if Initial_multiple_mapping in [0,1]:
            
            if sp.isLeaf:
                if sp.inital_ref==0:
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

                    sp.leftChild= driver(tr,sp.leftChild,sp_copy,sp_)
                    sp.rightChild =driver(tr,sp.rightChild,sp_copy,sp_)
                    #print('22',to_newick(sp))
                    
                    return sp
            else:
                    sp.evolve= 'Speciation'
                    sp.leftChild= driver(tr.leftChild,sp.leftChild,sp_copy,sp_)
                    sp.rightChild=driver(tr.rightChild,sp.rightChild,sp_copy,sp_)
                    
                    sp.cost=0
                    #print('1283',to_newick(sp))
                    
                    return sp
           



        else:
                

                #print(5)
                new_topo,cost =(ILS(tr_copy_1,sp_1,sp_1,Initial_multiple_mapping))
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


                recon.clean_up()
                recon.total_cost_()


                

               
                new_topo.label_internal()
                new_topo.map_recon(recon_1)
                recon_1.clean_up()
                recon_1.total_cost_()
                recon_1.cost= recon_1.cost+(Initial_multiple_mapping- cost)

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
                recon_1_cost =recon.cost
                recon_left.reset()
                recon_right.reset()
                #print(recon_1_cost)

                #print(recon_1.cost)
 
                if  recon_1.cost<recon_1_cost and  sp.isLeaf==None:
                        #print('NNI')
                        #print(to_newick(new_topo))
                        #print(to_newick(tr))
                        sp.refTo=[]
                        sp.evolve='NNI'
                        
                        if len(set(new_topo.leftChild.taxa).intersection(set(sp.leftChild.taxa)))>len(set(new_topo.rightChild.taxa).intersection(set(sp.leftChild.taxa))):
                           
                            sp.leftChild = driver(new_topo.leftChild,sp.leftChild,sp_copy,sp_) 
                            sp.rightChild = driver(new_topo.rightChild,sp.rightChild,sp_copy,sp_)
                        else:           
                            sp.leftChild = driver(new_topo.rightChild,sp.leftChild,sp_copy,sp_) 
                            sp.rightChild = driver(new_topo.leftChild,sp.rightChild,sp_copy,sp_)
                            
                        sp.cost=Initial_multiple_mapping- cost
                        #print('NNI',to_newick(sp))
                    
                        return sp

                if  recon_1.cost>=recon_1_cost or  sp.isLeaf:

                        #print('Duplication')
                        
                        sp.refTo=[]

                        recon_left = copy.deepcopy(sp)
                        recon_right = copy.deepcopy(sp)
                        recon_left.reset()
                        recon_right.reset()
                        recon_left.parent=sp
                        recon_right.parent=sp
                        sp.evolve='Duplication'
                        sp.refTo=sp.refTo[1:]

                        if sp.isLeaf:
                            #print(1)
                            if Initial_multiple_mapping==2:
                                recon_left.isLeaf=True
                                recon_right.isLeaf=True
                                recon_left.evolve='Speciation'

                                recon_right.evolve='Speciation'
                                sp.leftChild=recon_left
                                sp.rightChild=recon_right
                                return sp
                            #print(12455)
                            


                            tr.leftChild.reset()
                            tr.rightChild.reset()
                            recon_left.isLeaf=True
                            recon_right.isLeaf=True
                            clearcost(recon_left)
                            #print('11',to_newick(recon_left))

                            #print('11',to_newick(recon_right))

                            clearcost(recon_right)
                                                    

                            recon_left.optimize_cost(recon_left,tr.leftChild)
                            recon_right.optimize_cost(recon_right,tr.rightChild)

  
                                                          

                                                  
                            sp.leftChild = driver(tr.leftChild,recon_left,sp_copy,sp_) 
                            sp.rightChild = driver(tr.rightChild,recon_right,sp_copy,sp_)



 
                            return sp
                        if Initial_multiple_mapping==2:
                  
                            tr.leftChild.reset()
                            tr.rightChild.reset()

     
                            recon_left.optimize_cost(recon_left,tr.leftChild)
                            recon_right.optimize_cost(recon_right,tr.rightChild)

                            clearcost(recon_left)

                            clearcost(recon_right)
                            #print('12',to_newick(recon_left))

                            #print('12',to_newick(recon_right))
                                                  
                            sp.leftChild = driver(tr.leftChild,recon_left,sp_copy,sp_) 
                            sp.rightChild = driver(tr.rightChild,recon_right,sp_copy,sp_)

                            clearcost(sp)
                            #print('dups 22',to_newick(sp))
                            return sp

                        else:
                            tr.leftChild.reset()
                            tr.rightChild.reset()

                              
  
                                                        
   
                            recon_left.optimize_cost(recon_left,tr.leftChild)
                            recon_right.optimize_cost(recon_right,tr.rightChild)

                            #print('12',to_newick(recon_left))

                            #print('12',to_newick(recon_right))
                      
                            clearcost(recon_left)

                            clearcost(recon_right)

                            sp.leftChild = driver(tr.leftChild,recon_left,sp_copy,sp_) 
                            sp.rightChild = driver(tr.rightChild,recon_right,sp_copy,sp_)

                            clearcost(sp)
    
                            sp.cost=1
                            #print('dups 23',to_newick(sp))
                    
                            return sp
                else:

                        sp.evolve='Speciation'
                        sp.leftChild =driver(tr,sp.leftChild,sp_copy,sp_)
                        sp.rightChild = driver(tr,sp.rightChild,sp_copy,sp_)
                        return sp


