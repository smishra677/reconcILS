import sys
sys.path.append("../")
import utils.Tree as Tree
import copy
import utils.Tally as Tally
import argparse
import utils.ILS as ILS
import utils.readWrite as readWrite


a= Tree.Tree()
D_cost=  1
I_cost= 1
L_cost= 1


def clearid(sp,ori):
    import utils.idmaker_ as idmaker_
    if sp:

        sp.id =idmaker_.idmaker2().id
 

        clearid(sp.leftChild,ori)
            
        clearid(sp.rightChild,ori) 
            



            





def copy_event(re,sp):
    if sp:
        sp.NNI_.append(re.NNI_)
        copy_event(re.leftChild,sp.leftChild)
        copy_event(re.rightChild,sp.rightChild)





def setCost(sp):
    if sp:
        sp.inital_ref= len(sp.refTo)
    

        setCost(sp.leftChild)
        
        setCost(sp.rightChild)  




        
def reconcILS(tr,sp,sp_copy,sp_,visited):
    if sp:

        sp_copy = copy.deepcopy(sp)
        tr_copy_1 = copy.deepcopy(tr)

        sp_1 =copy.deepcopy(sp) 
        tr_copy_2 = copy.deepcopy(tr)
        

        Initial_multiple_mapping=len(sp.refTo)
        
        print('Multiple_mapping',Initial_multiple_mapping)




        if tr==None:
            if sp.isLeaf:
                if sp.inital_ref==0:
                    sp = sp.parse(sp.to_newick(sp))
                    sp.evolve='Loss'
                elif Initial_multiple_mapping>=2:
                    sp.evolve='Duplication'
                    sp = sp.parse(sp.to_newick(sp))
                return sp 
            else:
                if Initial_multiple_mapping>=2:
                    sp.evolve= 'Duplication'
                else:
                    sp.evolve='Speciation'
                return sp


        if Initial_multiple_mapping in [0,1]:

            if sp.isLeaf:
                if sp.inital_ref==0 and sp.parent.evolve!='Loss':
                    
                    sp = sp.parse(sp.to_newick(sp))
                    sp.evolve='Loss'
                    return sp
                else:
                    sp = sp.parse(sp.to_newick())
                    sp.evolve='Speciation'
                    return sp
            elif (sp.leftChild.isLeaf and sp.rightChild.isLeaf) :
                    if sp.evolve==None:
                        sp.evolve= 'Speciation'
                    if len(set(sp.leftChild.taxa).intersection(set(tr.taxa)))==0:
                        sp.leftChild.evolve='Loss'
                        sp.leftChild =sp.leftChild
                        sp.rightChild= reconcILS(tr,sp.rightChild,sp_copy,sp_,visited)
                        return sp
                    
                    if len(set(sp.rightChild.taxa).intersection(set(tr.taxa)))==0:
                        sp.rightChild.evolve='Loss'
                        sp.rightChild =sp.rightChild
                        sp.leftChild= reconcILS(tr,sp.leftChild,sp_copy,sp_,visited)
                        return sp
                    

                    if len(set(sp.leftChild.taxa).intersection(set(tr.leftChild.taxa)))>=len(set(sp.rightChild.taxa).intersection(set(tr.leftChild.taxa))):
                        sp.leftChild= reconcILS(tr.leftChild,sp.leftChild,sp_copy,sp_,visited)
                        sp.rightChild =reconcILS(tr.rightChild,sp.rightChild,sp_copy,sp_,visited)
                    else:
                        sp.leftChild= reconcILS(tr.rightChild,sp.leftChild,sp_copy,sp_,visited)
                        sp.rightChild =reconcILS(tr.leftChild,sp.rightChild,sp_copy,sp_,visited)                        
                    
                    
                    return sp
            
            elif sp.evolve!=None:
                if len(set(sp.leftChild.taxa).intersection(set(tr.leftChild.taxa)))>=len(set(sp.rightChild.taxa).intersection(set(tr.leftChild.taxa))):
                        sp.leftChild= reconcILS(tr.leftChild,sp.leftChild,sp_copy,sp_,visited)
                        sp.rightChild =reconcILS(tr.rightChild,sp.rightChild,sp_copy,sp_,visited)
                else:
                        sp.leftChild= reconcILS(tr.rightChild,sp.leftChild,sp_copy,sp_,visited)
                        sp.rightChild =reconcILS(tr.leftChild,sp.rightChild,sp_copy,sp_,visited)                          
                return sp
            else:
                    sp.evolve= 'Speciation'
                    if len(set(sp.leftChild.taxa).intersection(set(tr.taxa)))==0:
                        sp.leftChild.evolve='Loss'
                        sp.leftChild =sp.leftChild
                        sp.rightChild= reconcILS(tr,sp.rightChild,sp_copy,sp_,visited)
                        return sp
                    
                    if len(set(sp.rightChild.taxa).intersection(set(tr.taxa)))==0:
                        sp.rightChild.evolve='Loss'
                        sp.leftChild =reconcILS(tr,sp.leftChild,sp_copy,sp_,visited)
                        sp.rightChild= sp.rightChild
                        return sp
                    

                    if len(set(sp.leftChild.taxa).intersection(set(tr.leftChild.taxa)))>=len(set(sp.rightChild.taxa).intersection(set(tr.leftChild.taxa))):
                        sp.leftChild= reconcILS(tr.leftChild,sp.leftChild,sp_copy,sp_,visited)
                        sp.rightChild =reconcILS(tr.rightChild,sp.rightChild,sp_copy,sp_,visited)
                    else:
                        sp.leftChild= reconcILS(tr.rightChild,sp.leftChild,sp_copy,sp_,visited)
                        sp.rightChild =reconcILS(tr.leftChild,sp.rightChild,sp_copy,sp_,visited)                        
                                        
                    return sp
           



        else:
                if sp.isLeaf :
                    NNI_cost=1
                    duplication_cost=1
                else:

                    new_topo,cost =ILS.ILS().ILS(tr_copy_1,sp_1,sp_1,Initial_multiple_mapping,visited)

                    new_topo.reset()

                    recon_1 = copy.deepcopy(sp_1)
                    recon_1.reset()
                    recon_1.label_internal()


                    new_gene_tree =copy.deepcopy(tr_copy_2)
                    new_gene_tree.reset()
                    new_gene_tree.leftChild = sp.parse(new_gene_tree.leftChild.to_newick())
                    new_gene_tree.rightChild= sp.parse(new_gene_tree.rightChild.to_newick())


                    

                    recon_left = copy.deepcopy(sp)
                    recon_right = copy.deepcopy(sp)
                    
                    new_sp = copy.deepcopy(sp)
                    recon_right = sp.parse(recon_right.to_newick())
                    recon_left= sp.parse(recon_left.to_newick())


                    new_sp.leftChild=recon_left
                    new_sp.rightChild=recon_right

                    new_sp.reset()

                    recon_left.reset()
                    recon_right.reset()
                    recon_left.label_internal()
                    recon_right.label_internal()

                    
                    
                    
                    


                    new_gene_tree.leftChild.order_gene(recon_left)
                    new_gene_tree.rightChild.order_gene(recon_right)

                    new_gene_tree.leftChild.label_internal()

                    new_gene_tree.rightChild.label_internal()

                  
                    new_gene_tree.leftChild.map_gene(recon_left)

                    
                    new_gene_tree.rightChild.map_gene(recon_right)
                    recon_left.find_loss_sp(recon_left)
                    recon_right.find_loss_sp(recon_right)

                    


                    


                    recon_1_cost =recon_right.cost+recon_left.cost+1.1




                    recon_left.reset()

                    recon_right.reset()
                    
                    new_gene_tree.rightChild.reset()
                    new_gene_tree.leftChild.reset()
                    
                    new_topo.order_gene(recon_1)
                    new_topo.label_internal()
                    new_topo.map_gene(recon_1)
                    new_multiple= len(recon_1.refTo)

                    
                    recon_1.find_loss_sp(recon_1)
                    

                    #recon_1.cost= recon_1.cost+(Initial_multiple_mapping- cost)*1

                    NNI_cost= recon_1.cost+(Initial_multiple_mapping- cost)*1
                    duplication_cost=recon_1_cost



                if  NNI_cost<duplication_cost and  sp.isLeaf==None and (new_multiple<Initial_multiple_mapping or recon_1.cost==0):

                        sp.refTo=[]

                        new_topo.reset()
                        sp.clear_ref()


                        new_topo.order_gene(sp)

                                
                        new_topo.label_internal()
                    

                                        
                        new_topo.map_gene(sp)

                        
                        print(new_topo.to_newick())
                        

                        sp.cost=Initial_multiple_mapping- cost

                        if sp.evolve!=None:
                                if type(sp.evolve)==list:
                                    sp.evolve+=['NNI']
                                else:
                                    sp.evolve=[sp.evolve,'NNI']
                        else:
                            sp.evolve='NNI'
                        
                        
                        copy_event(sp_1,sp)

                        visited.append(new_topo.to_newick())
                        return reconcILS(new_topo,sp,sp_copy,sp_,visited)

                        
                       


                       

                if  NNI_cost>=duplication_cost or  sp.isLeaf or new_multiple>=Initial_multiple_mapping:

                        
                        recon_left = copy.deepcopy(sp)
                        clearid(recon_left,'Left')
                        
                        
                        
                        recon_right = copy.deepcopy(sp)
                        clearid(recon_right,'right')
                        
                        
                        recon_left.reset()
                        recon_right.reset()

                        recon_left.label_internal()
                        recon_right.label_internal()


                        sp.taxa=''



                       
                        if sp.evolve!=None:
                            if type(sp.evolve)==list:
                                sp.evolve+=['Duplication']
                            else:
                                sp.evolve=[sp.evolve,'Duplication']
                        else:
                            sp.evolve='Duplication'
                        

                        if 1==1:
                            tr.reset()


                            sp.clear_ref()

                            if tr.isLeaf:
                                tr.order_gene(recon_left)
                                tr.label_internal()
   
                                tr.map_gene(recon_left)
                                
                                sp.leftChild = reconcILS(tr,recon_left,sp_copy,sp_,visited) 
                                tr.reset()
                                tr.order_gene(recon_right)
                                tr.label_internal()
   
                                tr.map_gene(recon_right)
                                sp.rightChild = reconcILS(tr,recon_right,sp_copy,sp_,visited)
                                sp.children+=[sp.leftChild ]
                                sp.children+=[sp.rightChild]
                                sp.leftChild.parent=sp
                                sp.rightChild.parent=sp
                                
                                if sp.parent==None:
                                    sp.paralogy+=[(sp.id,sp.id)]
                                else:
                                    if sp.isLeaf:
                                        sp.paralogy+=[(sp.id, sp.id)]  
                                    else:  
                                        sp.paralogy+=[(sp.parent.id, sp.id)]  
                                    sp.paralogy+[(sp.parent.id, sp.id)]
                                sp.isLeaf=None
                            


                            
                            else:
                                tr.leftChild.order_gene(recon_left)
                                tr.rightChild.order_gene(recon_right)
                                tr.label_internal()

                                tr.leftChild.map_gene(recon_left)
                                sp.leftChild = reconcILS(tr.leftChild,recon_left,sp_copy,sp_,visited) 
                                tr.rightChild.map_gene(recon_right)
                                sp.rightChild = reconcILS(tr.rightChild,recon_right,sp_copy,sp_,visited)
                                sp.children+=[sp.leftChild ]
                                sp.children+=[sp.rightChild]
                                sp.leftChild.parent=sp
                                sp.rightChild.parent=sp

                                if sp.parent==None:
                                    sp.paralogy+=[(sp.id,sp.id)]
                                else:
                                    if sp.isLeaf:
                                        sp.paralogy+=[(sp.id, sp.id)]  
                                    else:  
                                        sp.paralogy+=[(sp.parent.id, sp.id)]
                                    
                                sp.isLeaf=None
                            
                           

                        
                            return sp
                        



def iterative_reconcILS(tr, sp, sp_copy, sp_, visited):
    stack = []
    eve=[]
    
    while True:
        if sp:

            sp_copy = copy.deepcopy(sp)
            tr_copy_1 = copy.deepcopy(tr)

            sp_1 = copy.deepcopy(sp)
            tr_copy_2 = copy.deepcopy(tr)

            Initial_multiple_mapping = len(sp.refTo)

            print('Multiple_mapping', Initial_multiple_mapping)

            if tr is None:
                if sp.isLeaf:
                    if sp.inital_ref == 0:
                        sp = sp.parse(sp.to_newick(sp))
                        sp.evolve = 'Loss'
                        
                    elif Initial_multiple_mapping >= 2:
                        sp.evolve = 'Duplication'
                        sp = sp.parse(sp.to_newick(sp))
                    eve.append(sp.evolve)
    
                else:
                    if Initial_multiple_mapping >= 2:
                        sp.evolve = 'Duplication'
                    else:
                        sp.evolve = 'Speciation'
                    eve.append(sp.evolve)
                    

            if Initial_multiple_mapping in [0,1]:

                if sp.isLeaf:
                    if sp.inital_ref==0 and sp.parent.evolve!='Loss':
                        
                        sp = sp.parse(sp.to_newick(sp))
                        sp.evolve='Loss'
                        eve.append(sp.evolve)
                        
                    else:
                        sp = sp.parse(sp.to_newick())
                        sp.evolve='Speciation'
                        eve.append(sp.evolve)
                        
                elif (sp.leftChild.isLeaf and sp.rightChild.isLeaf) :
                        if sp.evolve==None:
                            sp.evolve= 'Speciation'
                        if len(set(sp.leftChild.taxa).intersection(set(tr.taxa)))==0:
                            sp.leftChild.evolve='Loss'
                            eve.append(sp.leftChild.evolve )
                            
                            stack.append((tr,sp.rightChild))

                        
                        elif len(set(sp.rightChild.taxa).intersection(set(tr.taxa)))==0:
                            sp.rightChild.evolve='Loss'
                            eve.append(sp.rightChild.evolve )
                            stack.append((tr,sp.leftChild))
                        

                        elif len(set(sp.leftChild.taxa).intersection(set(tr.leftChild.taxa)))>=len(set(sp.rightChild.taxa).intersection(set(tr.leftChild.taxa))):
                            stack.append((tr.leftChild,sp.leftChild))
                            stack.append((tr.rightChild,sp.rightChild))
                            
                        else:
                            stack.append((tr.rightChild,sp.leftChild))
                            stack.append((tr.leftChild,sp.rightChild))
                            

                
                elif sp.evolve!=None:
                    if len(set(sp.leftChild.taxa).intersection(set(tr.leftChild.taxa)))>=len(set(sp.rightChild.taxa).intersection(set(tr.leftChild.taxa))):
                            stack.append((tr.leftChild,sp.leftChild))
                            stack.append((tr.rightChild,sp.rightChild))
                            
                    else:
                            stack.append((tr.rightChild,sp.leftChild))
                            stack.append((tr.leftChild,sp.rightChild))
                            
                else:
                        sp.evolve= 'Speciation'
                        if len(set(sp.leftChild.taxa).intersection(set(tr.taxa)))==0:
                        
                            sp.leftChild.evolve='Loss'
                            eve.append(sp.leftChild.evolve )
                            stack.append((tr,sp.rightChild))
                            
                        elif len(set(sp.rightChild.taxa).intersection(set(tr.taxa)))==0:
                            sp.rightChild.evolve='Loss'
                            eve.append(sp.rightChild.evolve )
                            stack.append((tr,sp.leftChild))
    

                        


                        elif len(set(sp.leftChild.taxa).intersection(set(tr.leftChild.taxa)))>=len(set(sp.rightChild.taxa).intersection(set(tr.leftChild.taxa))):
                            stack.append((tr.leftChild,sp.leftChild))
                            stack.append((tr.rightChild,sp.rightChild))
                            
                        else:
                            stack.append((tr.rightChild,sp.leftChild))
                            stack.append((tr.leftChild,sp.rightChild))
            



            else:
                    if sp.isLeaf :
                        NNI_cost=1
                        duplication_cost=1
                    else:

                        new_topo,cost =ILS.ILS().ILS(tr_copy_1,sp_1,sp_1,Initial_multiple_mapping,visited)
                        
                        new_topo.reset()

                        recon_1 = copy.deepcopy(sp_1)
                        recon_1.reset()
                        recon_1.label_internal()


                        new_gene_tree =copy.deepcopy(tr_copy_2)
                        new_gene_tree.reset()
                        new_gene_tree.leftChild = sp.parse(new_gene_tree.leftChild.to_newick())
                        new_gene_tree.rightChild= sp.parse(new_gene_tree.rightChild.to_newick())


                        

                        recon_left = copy.deepcopy(sp)
                        recon_right = copy.deepcopy(sp)
                        
                        new_sp = copy.deepcopy(sp)
                        recon_right = sp.parse(recon_right.to_newick())
                        recon_left= sp.parse(recon_left.to_newick())


                        new_sp.leftChild=recon_left
                        new_sp.rightChild=recon_right

                        new_sp.reset()

                        recon_left.reset()
                        recon_right.reset()
                        recon_left.label_internal()
                        recon_right.label_internal()

                        
                        
                        
                        


                        new_gene_tree.leftChild.order_gene(recon_left)
                        new_gene_tree.rightChild.order_gene(recon_right)

                        new_gene_tree.leftChild.label_internal()

                        new_gene_tree.rightChild.label_internal()

                    
                        new_gene_tree.leftChild.map_gene(recon_left)

                        
                        new_gene_tree.rightChild.map_gene(recon_right)
                        recon_left.find_loss_sp(recon_left)
                        recon_right.find_loss_sp(recon_right)

                        


                        


                        recon_1_cost =recon_right.cost+recon_left.cost+1.1




                        recon_left.reset()

                        recon_right.reset()
                        
                        new_gene_tree.rightChild.reset()
                        new_gene_tree.leftChild.reset()
                        
                        new_topo.order_gene(recon_1)
                        new_topo.label_internal()
                        new_topo.map_gene(recon_1)
                        new_multiple= len(recon_1.refTo)

                        
                        recon_1.find_loss_sp(recon_1)
                        

                        NNI_cost=recon_1.cost+(Initial_multiple_mapping- cost)*1

                        duplication_cost=recon_1_cost




                    if  NNI_cost<duplication_cost and  sp.isLeaf==None and (new_multiple<Initial_multiple_mapping or recon_1.cost==0) :
                            sp.refTo=[]
                            
                            new_topo.reset()
                            sp.clear_ref()


                            new_topo.order_gene(sp)

                                    
                            new_topo.label_internal()
                        

                                            
                            new_topo.map_gene(sp)

                            
                            
                            sp.cost=Initial_multiple_mapping- cost

         
                            copy_event(sp_1,sp)

                            visited.append(new_topo.to_newick())

                            eve.append(['NNI' for i in range(Initial_multiple_mapping- cost)])
                            stack.append((new_topo,sp))

                            
                        


                        

                    else:

                            
                            recon_left = copy.deepcopy(sp)
                            clearid(recon_left,'Left')
                            
                            
                            
                            recon_right = copy.deepcopy(sp)
                            clearid(recon_right,'right')
                            
                            
                            recon_left.reset()
                            recon_right.reset()

                            recon_left.label_internal()
                            recon_right.label_internal()


                            sp.taxa=''



                        
                        
                            eve.append('Duplication')
                            if 1==1:
                                tr.reset()


                                sp.clear_ref()

                                if tr.isLeaf:
                                    tr.order_gene(recon_left)
                                    tr.label_internal()
    
                                    tr.map_gene(recon_left)
                                    
                                    stack.append((tr,recon_left))
                                    tr.reset()
                                    tr.order_gene(recon_right)
                                    tr.label_internal()
    
                                    tr.map_gene(recon_right)
                                    stack.append((tr,recon_right))

                                    sp.leftChild=recon_left
                                    sp.rightChild=recon_right
                                    sp.children+=[sp.leftChild ]
                                    sp.children+=[sp.rightChild]
                                    
                                    sp.leftChild.parent=sp
                                    sp.rightChild.parent=sp
                                    
                                    if sp.parent==None:
                                        sp.paralogy+=[(sp.id,sp.id)]
                                    else:
                                        if sp.isLeaf:
                                            sp.paralogy+=[(sp.id, sp.id)]  
                                        else:  
                                            sp.paralogy+=[(sp.parent.id, sp.id)]  
                                        sp.paralogy+[(sp.parent.id, sp.id)]
                                    sp.isLeaf=None
                                


                                
                                else:
                                    tr.leftChild.order_gene(recon_left)
                                    tr.rightChild.order_gene(recon_right)
                                    tr.label_internal()

                                    tr.leftChild.map_gene(recon_left)
                                    stack.append((tr.leftChild,recon_left)) 
                                    tr.rightChild.map_gene(recon_right)
                                    stack.append((tr.rightChild,recon_right)) 
                                    
                                    sp.children+=[sp.leftChild ]
                                    sp.children+=[sp.rightChild]
                                

                                    if sp.parent==None:
                                        sp.paralogy+=[(sp.id,sp.id)]
                                    else:
                                        if sp.isLeaf:
                                            sp.paralogy+=[(sp.id, sp.id)]  
                                        else:  
                                            sp.paralogy+=[(sp.parent.id, sp.id)]
                                        
                                    sp.isLeaf=None
                                
                            




            if not stack:
                break
            else:
            
                
                tr, sp = stack.pop()

                

    return eve


def parse1():
    parser = argparse.ArgumentParser(description="reconcile Gene Tree with Species Tree")
    parser.add_argument('--spTree', type=str, help="Species_Tree")
    parser.add_argument('--gTree', type=str, help="gene_Tree")
    parser.add_argument('--output', type=str, help="Location and name to output")
    parser.add_argument('--D', type=str, help="Duplication Cost")
    parser.add_argument('--L', type=str, help="Loss Cost")
    parser.add_argument('--I', type=str, help="ILS Cost")
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
    
    if parser.D:
        D_cost=  parser.D
    if parser.I:
        I_cost=parser.I
    if parser.L:
        L_cost=parser.L
    
    red= readWrite.readWrite()
    tr= red.parse(gene_tree)

    tr=red.parse(red.to_newick(tr))


    sp=red.parse(sp_string)
    
    sp_copy= copy.deepcopy(sp)

    sp_copy_1= copy.deepcopy(sp)
    sp_copy.reset()

    tr.order_gene(sp)

    tr.label_internal()
    sp.label_internal()
    tr.map_gene(sp)


    #print('-----------------')




    setCost(sp)
    sp.isRoot=True
    tr.isRoot=True
    sp_copy.isRoot=True


    #reconcILS(tr,sp,sp_copy,sp,[])
    print(tr.to_newick())

    from multiprocessing import Process,Pool


    with Pool(4) as pool:
        eve = pool.starmap(iterative_reconcILS, [(tr,sp,sp_copy,sp,[])])
    #eve=iterative_reconcILS(tr,sp,sp_copy,sp,[])
    print(eve)

  
    #print('######################33')
    

    re_w = readWrite.readWrite()
    #li =re_w.sp_event(sp,[])
    print(sp.to_newick())
           
   
    #Tally.Tally().make_graph(sp,sp_string,gene_tree)
    



    sp.reset()

    tr.order_gene(sp)
    

    #print(Counter(li))
    dic={'Process':[],'Replicate':[],'Gene_tree':[],'Species_Tree':[],'Duplication':[],'NNI':[],'DLCILS':[],'Loss':[],'Hemiplasy':[],'RHemiplasy':[]}
    
    dic['Gene_tree']+=[tr.to_newick()]
    dic['Species_Tree']+=[sp_string]
        
    #print(li)
    dic= re_w.Create_pd('reconcILS',0,eve[0],dic)
    
    df = pd.DataFrame(dic)
    print(dic)
    exit()
    df.to_csv(parser.output, index=False)


if __name__ == "__main__":
    
    main()