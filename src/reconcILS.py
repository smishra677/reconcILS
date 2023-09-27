
import Tree
import copy
import Tally
import argparse
import ILS
import readWrite


a= Tree.Tree()
D_cost=  1
I_cost= 1
L_cost= 1


def clearid(sp,ori):
    import idmaker_
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
                
                #sp.leftChild= reconcILS(None,sp.leftChild,sp_copy,sp_)
                #sp.rightChild =reconcILS(None,sp.rightChild,sp_copy,sp_)
                return sp


        if Initial_multiple_mapping in [0,1]:

        
            if sp.isLeaf:
                if sp.inital_ref==0 and sp.parent.evolve!='Loss':
                    sp.evolve='Loss'

                    sp.cost=0
                    sp = sp.parse(sp.to_newick(sp))
                    sp.evolve='Loss'
                    return sp
                else:
                    sp.evolve='Speciation'
                    sp = sp.parse(sp.to_newick())
                    sp.evolve='Speciation'
                    return sp
            elif (sp.leftChild.isLeaf and sp.rightChild.isLeaf) :
                    if sp.evolve==None:
                        sp.evolve= 'Speciation'

                    if len(set(sp.leftChild.taxa).intersection(set(tr.taxa)))==0:
                        sp.leftChild.evolve='Loss'
                        print('right')
                        #clearid(sp,'right')
                        sp.leftChild =sp.leftChild
                        sp.rightChild= reconcILS(tr,sp.rightChild,sp_copy,sp_)
                        return sp
                    
                    if len(set(sp.rightChild.taxa).intersection(set(tr.taxa)))==0:
                        sp.rightChild.evolve='Loss'
                        print('left')
                        #clearid(sp,'Left')
                        sp.rightChild =sp.rightChild
                        sp.leftChild= reconcILS(tr,sp.leftChild,sp_copy,sp_)
                        return sp
                    

                    if len(set(sp.leftChild.taxa).intersection(set(tr.leftChild.taxa)))>=len(set(sp.rightChild.taxa).intersection(set(tr.leftChild.taxa))):


                        sp.leftChild= reconcILS(tr.leftChild,sp.leftChild,sp_copy,sp_)
                        sp.rightChild =reconcILS(tr.rightChild,sp.rightChild,sp_copy,sp_)
                    else:

                        sp.leftChild= reconcILS(tr.rightChild,sp.leftChild,sp_copy,sp_)
                        sp.rightChild =reconcILS(tr.leftChild,sp.rightChild,sp_copy,sp_)                        
                    
                    sp.cost=0
                    #print('1283',to_newick(sp))
                    
                    return sp
            elif sp.evolve!=None:
                if len(set(sp.leftChild.taxa).intersection(set(tr.leftChild.taxa)))>=len(set(sp.rightChild.taxa).intersection(set(tr.leftChild.taxa))):
                        sp.leftChild= reconcILS(tr.leftChild,sp.leftChild,sp_copy,sp_)
                        sp.rightChild =reconcILS(tr.rightChild,sp.rightChild,sp_copy,sp_)
                else:

                        sp.leftChild= reconcILS(tr.rightChild,sp.leftChild,sp_copy,sp_)
                        sp.rightChild =reconcILS(tr.leftChild,sp.rightChild,sp_copy,sp_)                        
                    
                    
                return sp
             
            else:
                    sp.evolve= 'Speciation'
      
                    if len(set(sp.leftChild.taxa).intersection(set(tr.taxa)))==0:
                        print('left1')
                        sp.leftChild.evolve='Loss'
                        #clearid(sp,'right')
                        sp.leftChild =sp.leftChild
                        sp.rightChild= reconcILS(tr,sp.rightChild,sp_copy,sp_)
                        return sp
                    
                    if len(set(sp.rightChild.taxa).intersection(set(tr.taxa)))==0:
                        print('right1')
                        sp.rightChild.evolve='Loss'
                        #clearid(sp,'Left')
                        sp.leftChild =reconcILS(tr,sp.leftChild,sp_copy,sp_)
                        sp.rightChild= sp.rightChild
                        return sp
                    

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
                if sp.isLeaf :
                    NNI_cost=1
                    duplication_cost=1
                else:

                    new_topo,cost =ILS.ILS().ILS(tr_copy_1,sp_1,sp_1,Initial_multiple_mapping)
                    #new_topo =parse(to_newick(new_topo))
                    new_topo.reset()

                    recon_1 = copy.deepcopy(sp_1)
                    recon_1.reset()
                    recon_1.label_internal()


                    new_gene_tree =copy.deepcopy(tr_copy_2)
                    new_gene_tree.reset()
                    new_gene_tree.leftChild = sp.parse(new_gene_tree.leftChild.to_newick())
                    new_gene_tree.rightChild= sp.parse(new_gene_tree.rightChild.to_newick())


                    

                    recon_left = copy.deepcopy(sp)
                    #recon_left.id =recon_left.id *10
                    recon_right = copy.deepcopy(sp)
                    #recon_right.id =recon_right.id *10
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

                    


                    
                    print(recon_left.cost)
                    print(recon_right.cost)


                    recon_1_cost =recon_right.cost+recon_left.cost+1
                    print(recon_1_cost)



                    recon_left.reset()

                    recon_right.reset()
                    
                    new_gene_tree.rightChild.reset()
                    new_gene_tree.leftChild.reset()
                    
                    new_topo.order_gene(recon_1)
                    new_topo.label_internal()
                    new_topo.map_gene(recon_1)

                    
                    recon_1.find_loss_sp(recon_1)
                    

                    recon_1.cost= recon_1.cost+(Initial_multiple_mapping- cost)*1
                    


                    NNI_cost=recon_1.cost
                    duplication_cost=recon_1_cost
    
 
                if  NNI_cost<=duplication_cost and  sp.isLeaf==None:
                        #print('NNI')
                        #print(to_newick(new_topo))
                        #print(to_newick(tr))
                        sp.refTo=[]
                        
                        new_topo.reset()
                        sp.clear_ref()


                        new_topo.order_gene(sp)

                                
                        new_topo.label_internal()
                    

                                        
                        new_topo.map_gene(sp)
                        print(sp.leftChild.refTo)
                        print(sp.rightChild.refTo)
                        #new_topo.id= new_topo.id*17
                        
                        sp.cost=Initial_multiple_mapping- cost
                        if sp.evolve!=None:
                            if type(sp.evolve)==list:
                                sp.evolve+=['NNI']
                            else:
                                sp.evolve=[sp.evolve,'NNI']
                        else:
                            sp.evolve='NNI'
                        #print('NNI',to_newick(new_topo))
                        copy_event(sp_1,sp)
                        #clearid(sp,'Left')
                        #clearid(sp,'Right')
                        

                        return reconcILS(new_topo,sp,sp_copy,sp_)

                        
                       


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

                if  NNI_cost>=duplication_cost or  sp.isLeaf:

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
                        if sp.evolve!=None:
                            if type(sp.evolve)==list:
                                sp.evolve+=['Duplication']
                            else:
                                sp.evolve=[sp.evolve,'Duplication']
                        else:
                            sp.evolve='Duplication'
                        #sp.refTo=sp.refTo[1:]
       
                        if 1==1:
                            print(1)

                            
    
                            



                            tr.reset()


                            sp.clear_ref()



                            
  

                            #clearcost_1(sp_copy,recon_left,recon_right)
                           
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
                                    sp.paralogy+=[(sp.id,sp.id)]
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
                                    sp.paralogy+=[(sp.id,sp.id)]
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


    print('-----------------')


    #
    print(red.to_newick(tr))

    setCost(sp)
    sp.isRoot=True
    tr.isRoot=True
    sp_copy.isRoot=True


    reconcILS(tr,sp,sp_copy,sp)
    #print('######################33')
    
    re_w = readWrite.readWrite()
    li =re_w.sp_event(sp,[])
           
   
    Tally.Tally().make_graph(sp,sp_string,gene_tree)
    



    sp.reset()

    tr.order_gene(sp)
    

    #print(Counter(li))
    dic={'Process':[],'Replicate':[],'Gene_tree':[],'Species_Tree':[],'Duplication':[],'NNI':[],'DLCILS':[],'Loss':[],'Hemiplasy':[],'RHemiplasy':[]}
    
    dic['Gene_tree']+=[tr.to_newick()]
    dic['Species_Tree']+=[sp_string]
        
    dic= re_w.Create_pd('reconcILS',0,li,dic)
    
    df = pd.DataFrame(dic)
    print(dic)
    df.to_csv(parser.output, index=False)


if __name__ == "__main__":
    
    main()