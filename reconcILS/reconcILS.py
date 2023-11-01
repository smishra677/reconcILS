import sys
sys.path.append("../")
import utils.Tree as Tree
import copy
import utils.Tally as Tally
import argparse
import utils.ILS as ILS
import utils.readWrite as readWrite


class reconcils:
    def __init__(self):
        self.a= Tree.Tree()
        self.D_cost=  1.1
        self.I_cost= 1
        self.L_cost= 1
        self.gene_tree=None
        self.address_dictionary=None
        self.V=False



    def clearid(self,sp,ori):
        import utils.idmaker_ as idmaker_
        if sp:

            sp.id =idmaker_.idmaker2().id
    

            self.clearid(sp.leftChild,ori)
                
            self.clearid(sp.rightChild,ori) 

    
    def  address_dict(self,sp_gene):
        if sp_gene:

            self.address_dict(sp_gene.leftChild)
            self.address_dict(sp_gene.rightChild)

         
         
     
                



                





    def copy_event(self,re,sp):
        if sp:
            #re.NNI_ = [lis for lis in re.NNI_ if lis != []]
            sp.NNI_.append(re.NNI_)
            #re.NNI_id_list = [lis for lis in re.NNI_id_list if lis != []]
            #sp.NNI_id_list.append(re.NNI_id_list)
            self.copy_event(re.leftChild,sp.leftChild)
            self.copy_event(re.rightChild,sp.rightChild)

    def search_event_(self,re,sp):
        if re:
            if sp.id==re.id:
                sp= re
                
                return
            self.search_event_(re.leftChild,sp)
            self.search_event_(re.rightChild,sp)



    def label_lost_child(self,tree):
        if tree:
            tree.event_list.append([-2,['L']])        
            self.label_lost_child(tree.leftChild)
            self.label_lost_child(tree.rightChild)


    def print_id(self,tree):
        if tree:
            print(tree.to_newick(),tree.id)
            self.print_id(tree.leftChild)
            self.print_id(tree.rightChild)
        


    def copy_loss(self,re,sp,new_topo):
        if sp:

            if sp.id==re.id:
                if re.parent==None:
                        co= Tree.Tree()
                        co.leftChild= copy.deepcopy(sp)
                        co.rightChild= copy.deepcopy(new_topo)

                       
                        sp.leftChild=co.leftChild
                        sp.rightChild=co.rightChild

                        co.leftChild.parent =sp
                        co.rightChild.parent =sp

                        sp.children=[sp.leftChild,sp.rightChild]
                        sp.isLeaf=None
            
                else:
                    if sp.parent.leftChild==sp:
                        co= Tree.Tree()
                        co.leftChild= copy.deepcopy(sp)
                        co.rightChild= copy.deepcopy(new_topo)
                        self.label_lost_child(co.rightChild)
                        

                        co.leftChild.parent =co
                        co.rightChild.parent =co

                        sp.parent.leftChild=co
                        sp.parent.isLeaf=None
                        
                        sp.parent.rightChild=sp.parent.rightChild
                        sp.parent.isLeaf=None
                        sp.parent.rightChild.parent=sp.parent
                        sp.parent.leftChild.parent=sp.parent
                       
                     

                        sp.event_list=[]
                        sp.parent.children=[sp.parent.leftChild,sp.parent.rightChild]
                        

                    else:
                        co= Tree.Tree()
                        co.leftChild= copy.deepcopy(sp)
                        co.rightChild= copy.deepcopy(new_topo)

                      
                        self.label_lost_child(co.rightChild)
                        
                        co.leftChild.parent =co
                        co.rightChild.parent =co

                        sp.parent.rightChild=co
                        sp.parent.leftChild=sp.parent.leftChild
                        sp.parent.isLeaf=None
                        sp.parent.rightChild.parent=sp.parent
                        sp.parent.leftChild.parent=sp.parent
                       
                     
                    
                        




                        sp.event_list=[]


                        sp.parent.children=[sp.parent.leftChild,sp.parent.rightChild]

                    sp.taxa=''

                    return
            else:
                self.copy_loss(re,sp.leftChild,new_topo)
                self.copy_loss(re,sp.rightChild,new_topo)           



    def copy_event_(self,re,sp,new_topo):
        if sp:
            if sp.id==re.id:
                    sp.event_list.append(re.event_list)
                    return
            else:
                self.copy_event_(re,sp.leftChild,new_topo)
                self.copy_event_(re,sp.rightChild,new_topo)
        
           



    def setCost(self,sp):
        if sp:
            sp.inital_ref= len(sp.refTo)
        

            self.setCost(sp.leftChild)
            
            self.setCost(sp.rightChild)  




    def reconcILS(self,tr,sp,sp_copy,sp_,visited):

        if sp:
            sp_copy = copy.deepcopy(sp)
            tr_copy_1 = copy.deepcopy(tr)

            sp_1 =copy.deepcopy(sp) 
            tr_copy_2 = copy.deepcopy(tr)
            

            Initial_multiple_mapping=len(sp.refTo)
            if self.V:
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
                            self.copy_loss(tr,self.gene_tree,sp.leftChild)
                            sp.rightChild= self.reconcILS(tr,sp.rightChild,sp_copy,sp_,visited)
                            return sp
                        
                        if len(set(sp.rightChild.taxa).intersection(set(tr.taxa)))==0:
                            sp.rightChild.evolve='Loss'
                            sp.rightChild =sp.rightChild
                            self.copy_loss(tr,self.gene_tree,sp.rightChild)
                            sp.leftChild= self.reconcILS(tr,sp.leftChild,sp_copy,sp_,visited)
                            return sp
                        

                        if len(set(sp.leftChild.taxa).intersection(set(tr.leftChild.taxa)))>=len(set(sp.rightChild.taxa).intersection(set(tr.leftChild.taxa))):
                            sp.leftChild= self.reconcILS(tr.leftChild,sp.leftChild,sp_copy,sp_,visited)
                            sp.rightChild =self.reconcILS(tr.rightChild,sp.rightChild,sp_copy,sp_,visited)
                        else:
                            sp.leftChild= self.reconcILS(tr.rightChild,sp.leftChild,sp_copy,sp_,visited)
                            sp.rightChild =self.reconcILS(tr.leftChild,sp.rightChild,sp_copy,sp_,visited)                        
                        

                        return sp
                
                elif sp.evolve!=None:
                    if len(set(sp.leftChild.taxa).intersection(set(tr.leftChild.taxa)))>=len(set(sp.rightChild.taxa).intersection(set(tr.leftChild.taxa))):
                            sp.leftChild= self.reconcILS(tr.leftChild,sp.leftChild,sp_copy,sp_,visited)
                            sp.rightChild =self.reconcILS(tr.rightChild,sp.rightChild,sp_copy,sp_,visited)
                    else:
                            sp.leftChild= self.reconcILS(tr.rightChild,sp.leftChild,sp_copy,sp_,visited)
                            sp.rightChild =self.reconcILS(tr.leftChild,sp.rightChild,sp_copy,sp_,visited)                          
                    return sp
                else:
                        sp.evolve= 'Speciation'
                        if len(set(sp.leftChild.taxa).intersection(set(tr.taxa)))==0:
                            sp.leftChild.evolve='Loss'
                            sp.leftChild =sp.leftChild
                            self.copy_loss(tr,self.gene_tree,sp.leftChild)
                            sp.rightChild= self.reconcILS(tr,sp.rightChild,sp_copy,sp_,visited)
                            return sp
                        
                        if len(set(sp.rightChild.taxa).intersection(set(tr.taxa)))==0:
                            sp.rightChild.evolve='Loss'
                            sp.rightChild= sp.rightChild
                            self.copy_loss(tr,self.gene_tree,sp.rightChild)
                            sp.leftChild =self.reconcILS(tr,sp.leftChild,sp_copy,sp_,visited)

                            
                            return sp
                        

                        if len(set(sp.leftChild.taxa).intersection(set(tr.leftChild.taxa)))>=len(set(sp.rightChild.taxa).intersection(set(tr.leftChild.taxa))):
                            sp.leftChild= self.reconcILS(tr.leftChild,sp.leftChild,sp_copy,sp_,visited)
                            sp.rightChild =self.reconcILS(tr.rightChild,sp.rightChild,sp_copy,sp_,visited)
                        else:
                            sp.leftChild= self.reconcILS(tr.rightChild,sp.leftChild,sp_copy,sp_,visited)
                            sp.rightChild =self.reconcILS(tr.leftChild,sp.rightChild,sp_copy,sp_,visited)                        
                    
                        return sp
            



            else:
                    if sp.isLeaf :
                        NNI_cost=1
                        duplication_cost=1
                    else:

                        new_topo,cost,bi_cos,child_ =ILS.ILS().ILS(tr_copy_1,sp_1,sp_1,Initial_multiple_mapping,[])
                    
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

                        


                        


                        recon_1_cost =recon_right.cost+recon_left.cost+self.D_cost




                        recon_left.reset()

                        recon_right.reset()
                        
                        new_gene_tree.rightChild.reset()
                        new_gene_tree.leftChild.reset()
                        
                        new_topo.order_gene(recon_1)
                        new_topo.label_internal()
                        new_topo.map_gene(recon_1)
                        new_multiple= len(recon_1.refTo)

                        
                        recon_1.find_loss_sp(recon_1)
                        

                        recon_1.cost= recon_1.cost+(Initial_multiple_mapping- cost)*self.I_cost

                        NNI_cost=recon_1.cost
                        duplication_cost=recon_1_cost


                

                    if  NNI_cost<duplication_cost and  sp.isLeaf==None and (new_multiple<Initial_multiple_mapping or bi_cos==0) :
                
                            sp.refTo=[]
                            
                            new_topo.reset()
                            sp.clear_ref()


                            new_topo.order_gene(sp)

                                    
                            new_topo.label_internal()
                        

                                            
                            new_topo.map_gene(sp)
                            
                            
                            
                            
                            sp.cost=1

                            if sp.evolve!=None:
                                if type(sp.evolve)==list:
                                    sp.evolve+=['NNI' for i in range(Initial_multiple_mapping- cost)]
                                    sp.cost=1
                                else:
                                    sp.evolve=[sp.evolve]
                                    sp.evolve+=['NNI' for i in range(Initial_multiple_mapping- cost)]
                                    sp.cost=1
                            else:
                                sp.evolve=['NNI' for i in range(Initial_multiple_mapping- cost)]

                            
                            


                            



                            
                            for i in list(child_):
                                if type(i[1].numbered_taxa)==set:
                                    li= '-'.join(list(i[1].numbered_taxa))
                                else:
                                    li =i[1].numbered_taxa
                                li.replace(';','')


                        
                                i[0].event_list+=[li,['I' for i in range(Initial_multiple_mapping- cost)]]

                                self.copy_event_(i[0],self.gene_tree,new_topo)

                        
                            


                            

                            

                            visited.append(new_topo.to_newick())
                            
                            return self.reconcILS(new_topo,sp,sp_copy,sp_,visited)

                            
                        


                        

                    if  NNI_cost>=duplication_cost or  sp.isLeaf or new_multiple>=Initial_multiple_mapping:

                            
                            recon_left = copy.deepcopy(sp)
                            
                            self.clearid(recon_left,'Left')
                            
                            
                            
                            recon_right = copy.deepcopy(sp)
                            self.clearid(recon_right,'right')
                            
                            
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


                            tr.event_list+=[-1,['D']]

                            self.copy_event_(tr,self.gene_tree,tr)




                           



                            
                            
        
                            if 1==1:
                                tr.reset()


                                sp.clear_ref()

                                if tr.isLeaf:
                                    tr.order_gene(recon_left)
                                    tr.label_internal()
    
                                    tr.map_gene(recon_left)
                                    
                                    sp.leftChild = self.reconcILS(tr,recon_left,sp_copy,sp_,visited) 
                                    tr.reset()
                                    tr.order_gene(recon_right)
                                    tr.label_internal()
    
                                    tr.map_gene(recon_right)
                                    sp.rightChild = self.reconcILS(tr,recon_right,sp_copy,sp_,visited)
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
                                    sp.leftChild = self.reconcILS(tr.leftChild,recon_left,sp_copy,sp_,visited) 
                                    tr.rightChild.map_gene(recon_right)
                                    sp.rightChild = self.reconcILS(tr.rightChild,recon_right,sp_copy,sp_,visited)
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
                            



    def iterative_reconcILS(self,tr, sp, sp_copy, sp_, visited):
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
                            eve.append('Loss')
                            
                        elif Initial_multiple_mapping >= 2:
                            sp.evolve = 'Duplication'
                            sp = sp.parse(sp.to_newick(sp))
                            eve.append('Duplication')

        
                    else:
                        if Initial_multiple_mapping >= 2:
                            sp.evolve = 'Duplication'
                            eve.append('Duplication')
                            
                        else:
                            sp.evolve = 'Speciation'
                            eve.append('Speciation')
                        

                if Initial_multiple_mapping in [0,1]:

                    if sp.isLeaf:
                        if sp.inital_ref==0 and sp.parent.evolve!='Loss':
                            
                            sp = sp.parse(sp.to_newick(sp))
                            sp.evolve='Loss'
                            eve.append('Loss')
                            
                        else:
                            sp = sp.parse(sp.to_newick())
                            sp.evolve='Speciation'
                            eve.append('Speciation')
                            
                    elif (sp.leftChild.isLeaf and sp.rightChild.isLeaf) :
                            if sp.evolve==None:
                                sp.evolve= 'Speciation'
                            if len(set(sp.leftChild.taxa).intersection(set(tr.taxa)))==0:
                                sp.leftChild.evolve='Loss'

                                eve.append('Loss')
                                stack.append((tr,sp.rightChild))

                            
                            elif len(set(sp.rightChild.taxa).intersection(set(tr.taxa)))==0:
                                sp.rightChild.evolve='Loss'
                                eve.append('Loss')
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
                                eve.append('Loss')
                            
                                stack.append((tr,sp.rightChild))
                                
                            elif len(set(sp.rightChild.taxa).intersection(set(tr.taxa)))==0:
                                sp.rightChild.evolve='Loss'
                                eve.append('Loss')
                        
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

                            new_topo,cost,bi_cos,child_=ILS.ILS().ILS(tr_copy_1,sp_1,sp_1,Initial_multiple_mapping,visited)
                            
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

                            


                            


                            recon_1_cost =recon_right.cost+recon_left.cost+self.D_cost




                            recon_left.reset()

                            recon_right.reset()
                            
                            new_gene_tree.rightChild.reset()
                            new_gene_tree.leftChild.reset()
                            
                            new_topo.order_gene(recon_1)
                            new_topo.label_internal()
                            new_topo.map_gene(recon_1)
                            new_multiple= len(recon_1.refTo)

                            
                            recon_1.find_loss_sp(recon_1)
                            

                            NNI_cost=recon_1.cost+(Initial_multiple_mapping- cost)*self.I_cost

                            duplication_cost=recon_1_cost


                            print('to_new',new_topo.to_newick())

                        if  NNI_cost<duplication_cost and  sp.isLeaf==None and (new_multiple<Initial_multiple_mapping or bi_cos==0):
                                sp.refTo=[]
                                
                                new_topo.reset()
                                sp.clear_ref()


                                new_topo.order_gene(sp)

                                        
                                new_topo.label_internal()
                            

                                                
                                new_topo.map_gene(sp)

                                
                                
                                sp.cost=Initial_multiple_mapping- cost

            
                                self.copy_event(sp_1,sp)

                                visited.append(new_topo.to_newick())

                                eve+=['NNI' for i in range(Initial_multiple_mapping- cost)]
                                stack.append((new_topo,sp))

                                
                            


                            

                        else:

                                
                                recon_left = copy.deepcopy(sp)
                                self.clearid(recon_left,'Left')
                                
                                
                                
                                recon_right = copy.deepcopy(sp)
                                self.clearid(recon_right,'right')
                                
                                
                                recon_left.reset()
                                recon_right.reset()

                                recon_left.label_internal()
                                recon_right.label_internal()


                                #sp.taxa=''



                            
                            
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
    parser.add_argument('--output', type=str, help="Location and name to output csv")
    parser.add_argument('--D', type=str, help="Duplication Cost")
    parser.add_argument('--L', type=str, help="Loss Cost")
    parser.add_argument('--I', type=str, help="ILS Cost")
    parser.add_argument('--V', type=str, help="Verbose Mode")
    args= parser.parse_args()
    return(args)




def main():
    reconcILS= reconcils()
    parser = parse1()
    from collections import Counter
    import pandas as pd
    import igraph as ig
    import matplotlib.pyplot as plt
    


    sp_string=parser.spTree
    gene_tree=parser.gTree
    
    if parser.D:
        reconcILS.D_cost=float(parser.D)
    if parser.I:
        reconcILS.I_cost=float(parser.I)
    if parser.L:
        reconcILS.L_cost=float(parser.L)
    if parser.V:
        if int(parser.V)==1:
            reconcILS.V=True

    
    red= readWrite.readWrite()
    tr= red.parse(gene_tree)
    sp=red.parse(sp_string)
    
    sp_copy= copy.deepcopy(sp)



    tr.order_gene(sp)
    tr.label_internal()
    sp.label_internal()
    tr.map_gene(sp)


    reconcILS.setCost(sp)
    sp.isRoot=True
    tr.isRoot=True
    sp_copy.isRoot=True
    

    reconcILS.gene_tree= copy.deepcopy(tr)


    reconcILS.reconcILS(tr,sp,sp_copy,sp,[])


    print(red.to_newick(reconcILS.gene_tree))


    #print('######################33')
    

    re_w = readWrite.readWrite()
 
    li =re_w.sp_event(sp,[])
    


    sp.reset()

    tr.order_gene(sp)
    

    #print(Counter(li))
    dic={'Process':[],'Replicate':[],'Gene_tree':[],'Species_Tree':[],'Duplication':[],'NNI':[],'DLCILS':[],'Loss':[],'Hemiplasy':[],'RHemiplasy':[]}
    
    dic['Gene_tree']+=[red.to_newick(reconcILS.gene_tree)]
    dic['Species_Tree']+=[sp_string]
        
    #print(li)
    dic= re_w.Create_pd('reconcILS',0,li,dic)
    
    df = pd.DataFrame(dic)
    
    
    df.to_csv(parser.output, index=False)


if __name__ == "__main__":
    
    main()