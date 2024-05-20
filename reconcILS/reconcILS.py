import sys
sys.path.append("../")
import utils.Tree as Tree
import copy
import utils.Tally as Tally
import argparse
import utils.ILS as ILS
import utils.readWrite as readWrite
import pickle
import time
import pandas as pd
import string
import json
from concurrent.futures import ThreadPoolExecutor, TimeoutError

sys.setrecursionlimit(50000)


class reconcils:
    def __init__(self):
        self.a= Tree.Tree()
        self.D_cost=  1.1
        self.I_cost= 1
        self.L_cost= 1
        self.gene_tree=None
        self.address_dictionary=None
        self.V=False
        self.F= False
        self.introgression=False
        self.time_ils=600
        self.thread=20

    def get_species(self,sp):
        flag1=False
        flag2=False
        length=[]
        return_it=[]
        stack = [sp]
        while stack:
                curr=stack.pop()

                
                

                if str(curr.taxa)!='0':
                    length.append(len(str(curr.taxa)))
                    return_it.append(str(curr.numbered_taxa))

                if curr.leftChild:
                    stack.append(curr.leftChild)
                if curr.rightChild:
                    stack.append(curr.rightChild)

        flag2  = all(len_ < 3 for len_ in length)
        flag1  = all(len_ > 3 for len_ in length)

        return return_it,not(flag1 or flag2),flag2


    def generate_letter_combinations(self):
        letters = list(string.ascii_uppercase)
        for first in letters:
            if first not in ['D','I','L']:
                yield first
        for first in letters:
            if first not in ['D','I','L']:
                for second in letters:
                    if second not in ['D','I','L']:
                        yield first + second



    def generate_code_species(self,mapping,tree):
        for species, letter in mapping.items():
            tree = tree.replace(","+species+")", " "+letter+")")
            tree = tree.replace("("+species+",", "("+letter+",")
        return tree

    def generate_code_gene(self,mapping,tree):
        for species, letter in mapping.items():
            tree = tree.replace(" "+species+",", " "+letter+",")
            tree = tree.replace(" "+species+")", " "+letter+")")
            tree = tree.replace("("+species+",", "("+letter+",")
            tree = tree.replace(","+species+")", ","+letter+")")
            tree = tree.replace(" "+species+"-", " "+letter+"-")
            tree = tree.replace("-"+species+"-", "-"+letter+"-")
            tree = tree.replace("-"+species+")", "-"+letter+")")
        return tree

    def generate_species_codes(self,mapping,tree):
        for species, letter in mapping.items():
            tree = tree.replace(species, letter)
        return tree

    
    def clearid(self,sp,ori):
        import utils.idmaker_ as idmaker_
        if sp:

            sp.id =idmaker_.idmaker2().id
    

            self.clearid(sp.leftChild,ori)
                
            self.clearid(sp.rightChild,ori) 

    
    
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
    
    def copy_numbered_taxa(self,sp,new_sp):
        if sp:
            new_sp.id= sp.id
            if sp.isLeaf==None:
                new_sp.event_list =sp.event_list
            else:
                new_sp.numbered_taxa= sp.numbered_taxa

            if new_sp.leftChild and new_sp.rightChild:

                if (sp.leftChild.to_newick()== new_sp.leftChild.to_newick()) and (sp.rightChild.to_newick()== new_sp.rightChild.to_newick()):

                    self.copy_numbered_taxa(sp.leftChild,new_sp.leftChild)
                    self.copy_numbered_taxa(sp.rightChild,new_sp.rightChild) 
                else:
                    self.copy_numbered_taxa(sp.rightChild,new_sp.leftChild)
                    self.copy_numbered_taxa(sp.leftChild,new_sp.rightChild) 


    def copy_loss(self,re,sp,new_topo):
        red= readWrite.readWrite()
        if sp:

            if sp.id==re.id:
                if re.parent==None:
                        co= Tree.Tree()
                        co.leftChild= sp
                        co.rightChild= new_topo
                        
                        co= co.parse(co.to_newick())
 
                        self.copy_numbered_taxa(sp,co.rightChild)
                        

                        self.label_lost_child(co.leftChild)
                        
                       
                        self.gene_tree=co
                        sp.id=co.id

                        co.leftChild.parent =sp
                        co.rightChild.parent =sp

                        sp.children=[sp.leftChild,sp.rightChild]
                        sp.isLeaf=None
            
                else:
                    if sp.parent.leftChild==sp:
                            co= Tree.Tree()
                            co.leftChild= sp
                            co.rightChild= new_topo


                            co= co.parse(co.to_newick())

                            self.copy_numbered_taxa(sp,co.rightChild)
                            


                            
                            self.label_lost_child(co.leftChild)


                            sp.parent.leftChild=co
                            sp.parent.isLeaf=None
                            
                            sp.parent.rightChild=sp.parent.rightChild
                            sp.parent.isLeaf=None
                            sp.parent.rightChild.parent=sp.parent
                            sp.parent.leftChild.parent=sp.parent
                        
                        

                            #sp.event_list=[]
                            sp.parent.children=[sp.parent.leftChild,sp.parent.rightChild]
                        
                            
                            

                    else:
                    
                            co= Tree.Tree()
                            co.leftChild= sp
                            co.rightChild= new_topo
   
                            co= co.parse(co.to_newick())


                            self.copy_numbered_taxa(sp,co.rightChild)
                            

                        
                            self.label_lost_child(co.leftChild)
                            
  
                            sp.parent.rightChild=co
                            sp.parent.leftChild=sp.parent.leftChild
                            sp.parent.isLeaf=None
                            sp.parent.rightChild.parent=sp.parent
                            sp.parent.leftChild.parent=sp.parent
                       
                        

                            
                            sp.parent.children=[sp.parent.leftChild,sp.parent.rightChild]
                       
                     
                    
                        




                        #sp.event_list=[]


                        

                    sp.taxa=''

                    return
            else:
                self.copy_loss(re,sp.leftChild,new_topo)
                self.copy_loss(re,sp.rightChild,new_topo)           

          


    
    def copy_event_(self, re, sp, new_topo):
        stack = [sp] 
        
        while stack:
            current_node = stack.pop()  
            
            if current_node:
                
                if current_node.id == re.id:
                    current_node.event_list.append(re.event_list)
                    return 
                
               
                stack.append(current_node.leftChild)
                stack.append(current_node.rightChild)



    def get_edges(self,sp):
        return_it=[(repr(list(sorted(sp.taxa))),' to ',repr(list(sorted(sp.taxa))))]
        stack = [sp]
        while stack:
            curr=stack.pop()

            if curr.leftChild:
                if curr.leftChild.isLeaf:
                    return_it.append((repr(list(sorted(curr.taxa))),' to ',repr(list(sorted({curr.leftChild.taxa})))))
                    stack.append(curr.leftChild)
                else:
                    return_it.append((repr(list(sorted(curr.taxa))),' to ',repr(list(sorted(curr.leftChild.taxa)))))
                    stack.append(curr.leftChild)
            if curr.rightChild:
                if curr.rightChild.isLeaf:
                    return_it.append((repr(list(sorted(curr.taxa))),' to ',repr(list(sorted({curr.rightChild.taxa})))))
                    stack.append(curr.rightChild)
                else:
                    return_it.append((repr(list(sorted(curr.taxa))),' to ',repr(list(sorted(curr.rightChild.taxa)))))
                    stack.append(curr.rightChild)

        return return_it         

    
    def edge_to_event(self,sp,dic,flag):
        stack = [sp]
        if len(sp.event_list)<2:
            if sp.isLeaf:
                sp.event_list+=[[dic[(repr(list(sorted({sp.taxa}))),' to ',repr(list(sorted({sp.taxa}))))],'Up']]
            else:
                sp.event_list+=[[dic[(repr(list(sorted(sp.taxa))),' to ',repr(list(sorted(sp.taxa))))],'Up']]
        else:
            if sp.isLeaf:
                sp.event_list[flag][0]['D']+=dic[(repr(list(sorted({sp.taxa}))),' to ',repr(list(sorted({sp.taxa}))))]['D']
                sp.event_list[flag][0]['I']+=dic[(repr(list(sorted({sp.taxa}))),' to ',repr(list(sorted({sp.taxa}))))]['I']
                sp.event_list[flag][0]['L']+=dic[(repr(list(sorted({sp.taxa}))),' to ',repr(list(sorted({sp.taxa}))))]['L']
            else:
                sp.event_list[flag][0]['D']+=dic[(repr(list(sorted(sp.taxa))),' to ',repr(list(sorted(sp.taxa))))]['D']
                sp.event_list[flag][0]['I']+=dic[(repr(list(sorted(sp.taxa))),' to ',repr(list(sorted(sp.taxa))))]['I']
                sp.event_list[flag][0]['L']+=dic[(repr(list(sorted(sp.taxa))),' to ',repr(list(sorted(sp.taxa))))]['L']
        while stack:
            curr=stack.pop()
            if curr.parent:
                if len(curr.event_list)<2:
                    if curr.isLeaf:
                        curr.event_list+= [[dic[(repr(list(sorted(curr.parent.taxa))),' to ',repr(list(sorted({curr.taxa}))))],'Up']]
                    else:
                        curr.event_list+= [[dic[(repr(list(sorted(curr.parent.taxa))),' to ',repr(list(sorted(curr.taxa))))],'Up']]
                else:
                    if curr.isLeaf:
                        curr.event_list[flag][0]['D']+=dic[(repr(list(sorted(curr.parent.taxa))),' to ',repr(list(sorted({curr.taxa}))))]['D']
                        curr.event_list[flag][0]['I']+=dic[(repr(list(sorted(curr.parent.taxa))),' to ',repr(list(sorted({curr.taxa}))))]['I']
                        curr.event_list[flag][0]['L']+=dic[(repr(list(sorted(curr.parent.taxa))),' to ',repr(list(sorted({curr.taxa}))))]['L']

                    else:
                        curr.event_list[flag][0]['D']+=dic[(repr(list(sorted(curr.parent.taxa))),' to ',repr(list(sorted(curr.taxa))))]['D']
                        curr.event_list[flag][0]['I']+=dic[(repr(list(sorted(curr.parent.taxa))),' to ',repr(list(sorted(curr.taxa))))]['I']
                        curr.event_list[flag][0]['L']+=dic[(repr(list(sorted(curr.parent.taxa))),' to ',repr(list(sorted(curr.taxa))))]['L']

            if curr.leftChild:
                stack.append(curr.leftChild)
            if curr.rightChild:
                stack.append(curr.rightChild)
        

    
    def setCost(self,sp):
        if sp:
            sp.inital_ref= len(sp.refTo)
        

            self.setCost(sp.leftChild)
            
            self.setCost(sp.rightChild)  




    def reconcILS(self,tr,sp,sp_copy,sp_,visited):

        if sp:
            #print(sp.to_newick(),tr.to_newick())
           
            #sp_copy = copy.deepcopy(sp)

            serialized_instance = pickle.dumps(sp)
            sp_copy = pickle.loads(serialized_instance)

            #tr_copy_1 = copy.deepcopy(tr)

            serialized_instance = pickle.dumps(tr)
            tr_copy_1 = pickle.loads(serialized_instance)


            #sp_1 =copy.deepcopy(sp) 

            serialized_instance = pickle.dumps(sp)
            sp_1 = pickle.loads(serialized_instance)


            #tr_copy_2 = copy.deepcopy(tr)

            serialized_instance = pickle.dumps(tr)
            tr_copy_2 = pickle.loads(serialized_instance)


            Initial_multiple_mapping=len(sp.refTo)
            if self.V:
                print('Multiple_mapping',Initial_multiple_mapping)

          
            if sp.isLeaf==None:
                left_=set(sp.leftChild.taxa)
                right_= set(sp.rightChild.taxa)

                '''
                if sp.leftChild.isLeaf!=None:
                    left_=set({sp.leftChild.taxa})

                if sp.rightChild.isLeaf!=None:
                    right_= set({sp.rightChild.taxa})
                '''
            tr_=set(tr.taxa)
            '''
            if tr.isLeaf:
                tr_=set({tr.taxa})
            '''

                        






            if Initial_multiple_mapping in [0,1]:
                if sp.isLeaf:
                    if sp.inital_ref==0 and sp.parent.evolve!='Loss':
                        #print('111')
                        sp.evolve='Loss'
                        return sp
                    else:
                        sp.evolve='Speciation'
                        return sp
                elif (sp.leftChild.isLeaf and sp.rightChild.isLeaf) :
                        #print('2')
                        if sp.evolve==None:
                            sp.evolve= 'Speciation'
                        if len(left_.intersection(tr_))==0:
                            sp.leftChild.evolve='Loss'
                            sp.leftChild =sp.leftChild
                            #self.copy_loss(tr,self.gene_tree,sp.leftChild)
                            sp.rightChild= self.reconcILS(tr,sp.rightChild,sp_copy,sp_,visited)
                            return sp
                        
                        if len(right_.intersection(tr_))==0:
                            sp.rightChild.evolve='Loss'
                            sp.rightChild =sp.rightChild
                            #self.copy_loss(tr,self.gene_tree,sp.rightChild)
                            sp.leftChild= self.reconcILS(tr,sp.leftChild,sp_copy,sp_,visited)
                            return sp
                        

                        if len(left_.intersection(tr_left))>=len(right_.intersection(tr_left)):
                            sp.leftChild= self.reconcILS(tr.leftChild,sp.leftChild,sp_copy,sp_,visited)
                            sp.rightChild =self.reconcILS(tr.rightChild,sp.rightChild,sp_copy,sp_,visited)
                        else:
                            sp.leftChild= self.reconcILS(tr.rightChild,sp.leftChild,sp_copy,sp_,visited)
                            sp.rightChild =self.reconcILS(tr.leftChild,sp.rightChild,sp_copy,sp_,visited)                        
                        

                        return sp
                
                elif sp.evolve!=None:
                    if len(left_.intersection(tr_left))>=len(right_.intersection(tr_left)):
                            sp.leftChild= self.reconcILS(tr.leftChild,sp.leftChild,sp_copy,sp_,visited)
                            sp.rightChild =self.reconcILS(tr.rightChild,sp.rightChild,sp_copy,sp_,visited)
                    else:
                            sp.leftChild= self.reconcILS(tr.rightChild,sp.leftChild,sp_copy,sp_,visited)
                            sp.rightChild =self.reconcILS(tr.leftChild,sp.rightChild,sp_copy,sp_,visited)                          
                    return sp
                else:   
                        #print(113232)
                        sp.evolve= 'Speciation'


                        if len(left_.intersection(tr_))==0:
                            sp.leftChild.evolve='Loss'
                            sp.leftChild =sp.leftChild
                            #self.copy_loss(tr,self.gene_tree,sp.leftChild)
                            sp.rightChild= self.reconcILS(tr,sp.rightChild,sp_copy,sp_,visited)
                            return sp

                        #print('-----------------------------------------------')
                        
                        if len(right_.intersection(tr_))==0:
                            sp.rightChild.evolve='Loss'
                            sp.rightChild= sp.rightChild
                            #self.copy_loss(tr,self.gene_tree,sp.rightChild)
                            sp.leftChild =self.reconcILS(tr,sp.leftChild,sp_copy,sp_,visited)

                            
                            return sp
                        

                        if len(left_.intersection(tr_left))>=len(right_.intersection(tr_left)):
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

                        #recon_1 = copy.deepcopy(sp_1)

                        serialized_instance = pickle.dumps(sp_1)
                        recon_1 = pickle.loads(serialized_instance)

                        recon_1.reset()
                        recon_1.label_internal()


                        #new_gene_tree =copy.deepcopy(tr_copy_2)
                        serialized_instance = pickle.dumps(tr_copy_2)
                        new_gene_tree = pickle.loads(serialized_instance)


                        new_gene_tree.reset()
                        #new_gene_tree.leftChild = sp.parse(new_gene_tree.leftChild.to_newick())
                        #new_gene_tree.rightChild= sp.parse(new_gene_tree.rightChild.to_newick())


                        

                        #recon_left = copy.deepcopy(sp)
                        serialized_instance = pickle.dumps(sp)
                        recon_left = pickle.loads(serialized_instance)
                        
                        #recon_right = copy.deepcopy(sp)
                        serialized_instance = pickle.dumps(sp)
                        recon_right = pickle.loads(serialized_instance)
                        
                        
                        #new_sp = copy.deepcopy(sp)
                        serialized_instance = pickle.dumps(sp)
                        new_sp = pickle.loads(serialized_instance)


                        #recon_right = sp.parse(recon_right.to_newick())
                        #recon_left= sp.parse(recon_left.to_newick())


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

                        


                        


                        recon_1_cost =(recon_right.cost+recon_left.cost)*self.L_cost+self.D_cost
  
                        



                        recon_left.reset()

                        recon_right.reset()
                        
                        new_gene_tree.rightChild.reset()
                        new_gene_tree.leftChild.reset()
                        
                        new_topo.order_gene(recon_1)
                        new_topo.label_internal()
                        new_topo.map_gene(recon_1)
                        new_multiple= len(recon_1.refTo)

                        
                        recon_1.find_loss_sp(recon_1)
                        

                        recon_1.cost= recon_1.cost*self.L_cost+(Initial_multiple_mapping- cost)*self.I_cost

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


                        
                                #i[0].event_list+=[li,['I' for i in range(Initial_multiple_mapping- cost)]]
                                i[0].event_list+=[li,['I']]

                                #self.copy_event_(i[0],self.gene_tree,new_topo)

   

                            #visited.append(new_topo.to_newick())
                            
                            return self.reconcILS(new_topo,sp,sp_copy,sp_,visited)


                    if  NNI_cost>=duplication_cost or  sp.isLeaf or new_multiple>=Initial_multiple_mapping:
                            #print('Dups')
                            #recon_left = copy.deepcopy(sp)
                            serialized_instance = pickle.dumps(sp)
                            recon_left = pickle.loads(serialized_instance)

                            
                            self.clearid(recon_left,'Left')
                            
                            
                            
                            #recon_right = copy.deepcopy(sp)
                            serialized_instance = pickle.dumps(sp)
                            recon_right = pickle.loads(serialized_instance)

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

                            #self.copy_event_(tr,self.gene_tree,tr)
        
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


                                    sp.isLeaf=None
                                
                            

                            
                                return sp
                            


    
    def iterative_reconcILS(self,tr, sp, sp_copy, sp_, visited):
        stack = []
        eve=[]

        
        while True:
            if sp:

                Initial_multiple_mapping = len(sp.refTo)
                if self.V:
                    print('Multiple_mapping', Initial_multiple_mapping)
                    print('species',sp.to_newick())
                    print('gene',tr.to_newick())



                if sp.isLeaf==None:
                    left_=set(sp.leftChild.taxa)
                    right_= set(sp.rightChild.taxa)
                    
                    if sp.leftChild.isLeaf!=None:
                        left_=set({sp.leftChild.taxa})

                    if sp.rightChild.isLeaf!=None:
                        right_= set({sp.rightChild.taxa})
                    
                
                tr_=set(tr.taxa)
                
                if tr.isLeaf:
                    tr_=set({tr.taxa})
                
                tr_left=None
                if tr.leftChild:
                    tr_left=set(tr.leftChild.taxa)
                                
                    if tr.leftChild.isLeaf:
                        tr_left=set({tr.leftChild.taxa})


                if Initial_multiple_mapping in [0,1]:

                    if sp.isLeaf:
                        #print(1)
                        if sp.inital_ref==0 and sp.parent.evolve!='Loss':
                            #sp = sp.parse(sp.to_newick(sp))
                            sp.evolve='Loss'
                            eve.append('Loss')
                            
                        else:
                            #sp = sp.parse(sp.to_newick())
                            sp.evolve='Speciation'
                            eve.append('Speciation')
                            
                    elif (sp.leftChild.isLeaf and sp.rightChild.isLeaf) :
                            if sp.evolve==None:
                                sp.evolve= 'Speciation'
                            if len(left_.intersection(tr_))==0:
                                #print('loss')
                                sp.leftChild.evolve='Loss'
                                self.copy_loss(tr,self.gene_tree,sp.leftChild)
                                eve.append([sp.taxa, sp.leftChild.taxa,'L',sp.leftChild.isLeaf])
                                stack.append((tr,sp.rightChild))

                            
                            elif len(right_.intersection(tr_))==0:
                                sp.rightChild.evolve='Loss'
                                eve.append([sp.taxa, sp.rightChild.taxa,'L',sp.rightChild.isLeaf])
                                #print('loss1')
                                self.copy_loss(tr,self.gene_tree,sp.rightChild)
                                stack.append((tr,sp.leftChild))
                            

                            elif len(left_.intersection(tr_left))>=len(right_.intersection(tr_left)):
                                stack.append((tr.leftChild,sp.leftChild))
                                stack.append((tr.rightChild,sp.rightChild))
                                
                            else:
                                stack.append((tr.rightChild,sp.leftChild))
                                stack.append((tr.leftChild,sp.rightChild))
                                

                    
                    elif sp.evolve!=None:
                        if len(left_.intersection(tr_left))>=len(right_.intersection(tr_left)):
                                stack.append((tr.leftChild,sp.leftChild))
                                stack.append((tr.rightChild,sp.rightChild))
                                
                        else:
                                stack.append((tr.rightChild,sp.leftChild))
                                stack.append((tr.leftChild,sp.rightChild))
                                
                    else:
                            sp.evolve= 'Speciation'
                            if len(left_.intersection(tr_))==0:
                                #print('loss2')
                                self.copy_loss(tr,self.gene_tree,sp.leftChild)
                                sp.leftChild.evolve='Loss'
                                eve.append([sp.taxa, sp.leftChild.taxa,'L',sp.leftChild.isLeaf])
                            
                                stack.append((tr,sp.rightChild))
                                
                            elif len(right_.intersection(tr_))==0:
                                sp.rightChild.evolve='Loss'
                                eve.append([sp.taxa, sp.rightChild.taxa,'L',sp.rightChild.isLeaf])
                                #print('loss3')
                                self.copy_loss(tr,self.gene_tree,sp.rightChild)
                                #print('done')
                                stack.append((tr,sp.leftChild))
        

                            


                            elif len(left_.intersection(tr_left))>=len(right_.intersection(tr_left)):
                                stack.append((tr.leftChild,sp.leftChild))
                                stack.append((tr.rightChild,sp.rightChild))
                                
                            else:
                                stack.append((tr.rightChild,sp.leftChild))
                                stack.append((tr.leftChild,sp.rightChild))
                



                else:   
                        if sp.isLeaf :
                            NNI_cost=1
                            duplication_cost=1
                            recon_left,recon_right=0,0
                        else:
                            tr_copy_1,tr_copy_2 = tr.deepcopy_double()
                            #print(tr_copy_1.to_newick(),sp.to_newick(),Initial_multiple_mapping)
                            
    
                            sp_1 = sp.deepcopy_single()
                            
                            def call_ils_function():
                                return ILS.ILS().ILS(tr_copy_1, sp_1, sp_1, Initial_multiple_mapping, [])

                            num_threads = self.thread

                            timeout_duration = self.time_ils
                            with ThreadPoolExecutor(max_workers=num_threads) as executor:
                                future = executor.submit(call_ils_function)

                                try:

                                    new_topo, cost, bi_cos, child_ = future.result(timeout=timeout_duration)

                                except TimeoutError:
                                    print(f"The function call timed out after {timeout_duration} seconds.")
                                    return []

                            new_topo.reset()

                            recon_1 = sp_1.deepcopy_single()
                            recon_1.reset()
                            recon_1.label_internal()


                            new_gene_tree =tr_copy_2.deepcopy_single()
                            new_gene_tree.reset()
                            

                            

                            recon_left,recon_right = sp.deepcopy_double()
                            
                            



                            recon_left.reset()
                            recon_right.reset()
                            recon_left.label_internal()
                            recon_right.label_internal()

                            
                            
                            
                            

                            #print(new_gene_tree.to_newick(),sp.to_newick(),Initial_multiple_mapping)
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


                            #print('to_new',new_topo.to_newick())

                        if  NNI_cost<duplication_cost and  sp.isLeaf==None and (new_multiple<Initial_multiple_mapping or bi_cos==0):
                                sp.refTo=[]
                                
                                new_topo.reset()
                                sp.clear_ref()


                                new_topo.order_gene(sp)

                                        
                                new_topo.label_internal()
                            

                                                
                                new_topo.map_gene(sp)

                                
                                
                                sp.cost=Initial_multiple_mapping- cost

            
                                #self.copy_event(sp_1,sp)

                                #visited.append(new_topo.to_newick())

                                #print((Initial_multiple_mapping- cost))
                                #print(child_)


                                
                                stack.append((new_topo,sp))

                                if self.introgression:

                                    df1 = pd.read_csv('./discordance.csv', converters={'events': eval}) 

                                for i in list(child_):
                                    if type(i[1].numbered_taxa)==set:
                                        li= '-'.join(list(i[1].numbered_taxa))
                                    else:
                                        li =i[1].numbered_taxa

                                    li.replace(';','')



                            
                                    #i[0].event_list+=[li,['I' for i in range(Initial_multiple_mapping- cost)]]
                                    i[0].event_list+=[li,['I']]
                                    self.copy_event_(i[0],self.gene_tree,new_topo)
                                    i[1].label_internal()

                                    left_sp_t=sp.leftChild.taxa
                                    right_sp_t=sp.rightChild.taxa

                                    if sp.leftChild.isLeaf:
                                        left_sp_t={sp.leftChild.taxa}
                                    if sp.rightChild.isLeaf:
                                        right_sp_t={sp.rightChild.taxa}

                                    
                                    i_taxa= i[1].taxa
                                    
                                    if i[1].isLeaf:
                                        i_taxa={i[1].taxa}



                                            
                                    if len(i_taxa.intersection(left_sp_t))>=len(i_taxa.intersection(right_sp_t)):
                                        eve+=[[sp.taxa,sp.leftChild.taxa,'I',sp.leftChild.isLeaf]]
                                        if sp.leftChild.isLeaf:
                                            branch_name = [repr(list(sorted(sp.taxa))),' to ',repr(list(sorted({sp.leftChild.taxa})))]
                                        else:
                                            branch_name = [repr(list(sorted(sp.taxa))),' to ',repr(list(sorted(sp.leftChild.taxa)))]
                                    else:
                                        eve+=[[sp.taxa,sp.rightChild.taxa,'I',sp.rightChild.isLeaf]]
                                        if sp.rightChild.isLeaf:

                                            branch_name = [repr(list(sorted(sp.taxa))),' to ',repr(list(sorted({sp.rightChild.taxa})))]
                                        else:
                                            branch_name = [repr(list(sorted(sp.taxa))),' to ',repr(list(sorted(sp.rightChild.taxa)))]
                                        
                                            
                                        #self.copy_event_(i[0],self.gene_tree,new_topo) 
                                

                                if self.introgression:
                                    introgression =sorted([sorted(list(i[0].parent.leftChild.taxa)),sorted(list(i[0].parent.rightChild.taxa))])

                                    new_events = [''.join(introgression[0])+'-'+''.join(introgression[1])]       

        

                                    mask = df1['branch'] == repr(branch_name)
                                    red =readWrite.readWrite()

                                    df1.loc[mask, 'events'] = df1.loc[mask, 'events'].apply(lambda x: red.update_events(x, new_events))

                                    
                                                
                                    df1.to_csv('./discordance.csv', index=False)                               


                       

                        else:
                                if recon_left==0:
                                
                                    recon_left,recon_right = sp.deepcopy_double()
 
                                
                                self.clearid(recon_left,'Left')

                                self.clearid(recon_right,'right')
                                
                                
                                recon_left.reset()
                                recon_right.reset()

                                recon_left.label_internal()
                                recon_right.label_internal()


                                #sp.taxa=''



                                tr.event_list+=[-1,['D']]
                                self.copy_event_(tr,self.gene_tree,tr)
                                if sp.isRoot:
                                    eve.append([sp.taxa,sp.taxa,'D',sp.isLeaf])
                                else:

                                    parent= sp.parent
                                    while sp.taxa== parent.taxa:
                                        parent=parent.parent
                                        

                                    eve.append([parent.taxa,sp.taxa,'D',sp.isLeaf])

                                


                            
                            
                                
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
                                        sp.children=[sp.leftChild ]
                                        sp.children+=[sp.rightChild]
                                
                                        sp.leftChild.parent=sp
                                        sp.rightChild.parent=sp   
                                        sp.isLeaf=None
                                                 
                                    else:
                                        tr.leftChild.order_gene(recon_left)
                                        tr.rightChild.order_gene(recon_right)
                                        tr.label_internal()
                                        tr.leftChild.map_gene(recon_left)
                                        stack.append((tr.leftChild,recon_left)) 
                                        tr.rightChild.map_gene(recon_right)
                                        stack.append((tr.rightChild,recon_right)) 
                                        sp.leftChild=recon_left
                                        sp.rightChild=recon_right 
                                        sp.children=[sp.leftChild ]
                                        sp.children+=[sp.rightChild]
                                        sp.leftChild.parent=sp
                                        sp.rightChild.parent=sp
                                        sp.isLeaf=None
                                    
                                       
                if not stack:
                    
                    break
                else:
                    
                    
                    tr, sp = stack.pop()

                    

        return eve
    

    def read_trees(self,file):
        gene_tre= open(file)
        tr =gene_tre.read().strip().split('\n')
        gene_tre.close()
        return str(tr[0])



def parse1():
    parser = argparse.ArgumentParser(description="reconcile Gene Tree with Species Tree")
    parser.add_argument('--spTree', type=str, help="Species_Tree")
    parser.add_argument('--gTree', type=str, help="gene_Tree")
    parser.add_argument('--output', type=str, help="Location and name to output csv")
    parser.add_argument('--D', type=str, help="Duplication Cost")
    parser.add_argument('--L', type=str, help="Loss Cost")
    parser.add_argument('--I', type=str, help="ILS Cost")
    parser.add_argument('--V', type=str, help="Verbose Mode")
    parser.add_argument('--F', type=str, help="Read File Mode")
    parser.add_argument('--Delta', type=str, help="Record number of discordant topologies")
    args= parser.parse_args()
    return(args)




def main():
    reconcILS= reconcils()
    parser = parse1()
    from collections import Counter
    import pandas as pd
    import igraph as ig
    import matplotlib.pyplot as plt
    from concurrent.futures import ThreadPoolExecutor

    

    if parser.F:
        if int(parser.F)==1:
            reconcILS.F=True
    if parser.D:
        reconcILS.D_cost=float(parser.D)
    if parser.I:
        reconcILS.I_cost=float(parser.I)
    if parser.L:
        reconcILS.L_cost=float(parser.L)
    if parser.V:
        if int(parser.V)==1:
            reconcILS.V=True

    if parser.F:
        sp_string=reconcILS.read_trees(parser.spTree)
        gene_tree=reconcILS.read_trees(parser.gTree)
    else:
        sp_string=parser.spTree
        gene_tree=parser.gTree

    
    red= readWrite.readWrite()
    sp_string = sp_string.replace('e-', '0')
    gene_tree = gene_tree.replace('e-', '0')



    tr= red.parse_bio(gene_tree)
    sp=red.parse_bio(sp_string)

    species_names,flag1,flag2=reconcILS.get_species(sp)


    if flag1:
        print('Error: please make sure all species are named with more than 2> or less than 2<= characters')
        exit()

    
 

    if not flag2:
        map_ = reconcILS.generate_letter_combinations()


        mapping = {name: next(map_) for name in species_names}


        resolve_sp_string=reconcILS.generate_species_codes(mapping,red.to_newick(sp))

        resolved_gene_tree=reconcILS.generate_species_codes(mapping,red.to_newick(tr))

        inverse_mapping = {value: key for key, value in mapping.items()}



        file_path_name_dict = './'+parser.output[:-4]+'species_acr.json'
        file_path_inverse_name_dict = './'+parser.output[:-4]+'acr_species.json'

        with open(file_path_name_dict, 'w') as file:
            json.dump(mapping, file)

        with open(file_path_inverse_name_dict, 'w') as file:
            json.dump(inverse_mapping, file)
    else:
        resolved_gene_tree=red.to_newick(tr)
        resolve_sp_string =red.to_newick(sp)


    sp=red.parse(resolve_sp_string)
    tr=red.parse(resolved_gene_tree)


    
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
    re_w = readWrite.readWrite()
    species_edge_list=reconcILS.get_edges(sp)
    
    if reconcILS.introgression:
        red.write_introgression(sp)
    

    li= reconcILS.iterative_reconcILS(tr,sp,sp_copy,sp,[])
        
        
    
    #print('##############################################################################################')
    if len(li)==0:
        print('Error')
    else:   
        #print('##############################################################################################')
        li.reverse()
        dic_reconcILS={}
        for i in li:
            if type(i)==list:
                    if (repr(list(sorted(i[0]))),' to ',repr(list(sorted(i[1]))))  in dic_reconcILS.keys():
                        if i[3]==True:
                            dic_reconcILS[(repr(list(sorted(i[0]))),' to ',repr(list(sorted({i[1]}))))]+=[i[2]]
                        else:
                            dic_reconcILS[(repr(list(sorted(i[0]))),' to ',repr(list(sorted(i[1]))))]+=[i[2]]
                    else:   
                        if i[3]==True:
                            dic_reconcILS[(repr(list(sorted(i[0]))),' to ',repr(list(sorted({i[1]}))))]=[i[2]]
                        else:
                            dic_reconcILS[(repr(list(sorted(i[0]))),' to ',repr(list(sorted(i[1]))))]=[i[2]]

    #print('##############################################################################################')

    eve_rec= {'D':0,'I':0,'L':0}
    for i in species_edge_list:
            rec_in_dic= {'D':0,'I':0,'L':0}
            if i in dic_reconcILS:
                dic_in= dict(Counter(dic_reconcILS[i]))
                for j in rec_in_dic:
                        if j in dic_in:
                            rec_in_dic[j]=dic_in[j]
                            eve_rec[j]+=dic_in[j]
                dic_reconcILS[i]=rec_in_dic
            else:
                dic_reconcILS[i]={'D':0,'I':0,'L':0}
                                
  



                                
                                
    sp_small=red.parse(resolve_sp_string)
    sp_small.label_internal()


    reconcILS.edge_to_event(sp_small,dic_reconcILS,0) 

    dic={'Process':[],'Replicate':[],'Gene_tree':[],'Species_Tree':[],'Duplication':[],'NNI':[],'Loss':[]}




    if not flag2:
        dic['Gene_tree']+=[reconcILS.generate_code_gene(inverse_mapping,red.to_newick(reconcILS.gene_tree))]
        dic['Species_Tree']+=[reconcILS.generate_code_gene(inverse_mapping,red.to_newick_sp(sp_small))]

    else:
        dic['Gene_tree']+=[red.to_newick(reconcILS.gene_tree)]
        dic['Species_Tree']+=[red.to_newick_sp(sp_small)]


    dic['Replicate']+=[0]
    dic['Process']+=['reconcILS']  
    dic['Duplication']+=[eve_rec['D']]
    dic['NNI']+=[eve_rec['I']]
    dic['Loss']+=[eve_rec['L']]

    df = pd.DataFrame(dic)
    print(dic)

    df.to_csv(parser.output, index=False)

    dic_log ={'Gene_Tree':[tr.to_newick()],'Species_Tree':[sp_string],'Duplication_cost':[str(reconcILS.D_cost)],'NNI_cost':[str(reconcILS.I_cost)],'Loss_cost':[str(reconcILS.L_cost)]}
    df_log=pd.DataFrame(dic_log)
    df_log.to_csv(parser.output[:-4]+'_log.csv',index=False)
if __name__ == "__main__":
    
    main()
