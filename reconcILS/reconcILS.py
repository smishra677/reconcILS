import sys
sys.path.append("../")
import utils.Tree as Tree
import copy
import utils.Tally as Tally
import argparse
import utils.ILS as ILS
import utils.readWrite as readWrite
import pickle
import gc

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
        


    def copy_loss(self, re, sp, new_topo):
        stack = [(re, sp)]
        while stack:
            current_re, current_sp = stack.pop()

            if not current_sp:
                continue

            if current_sp.id == current_re.id:
                print(1)
                if current_re.parent is None:
                    co = Tree.Tree()
                    print(2)

                    serialized_instance = pickle.dumps(current_sp)
                    co.leftChild = pickle.loads(serialized_instance)

                    serialized_instance = pickle.dumps(new_topo)
                    co.rightChild = pickle.loads(serialized_instance)

                    self.label_lost_child(co.rightChild)

                    current_sp.leftChild = co.leftChild
                    current_sp.rightChild = co.rightChild
                    current_sp.id = co.id

                    co.leftChild.parent = current_sp
                    co.rightChild.parent = current_sp

                    current_sp.children = [current_sp.leftChild, current_sp.rightChild]
                    current_sp.isLeaf = None
                else:
                    if current_sp.parent.leftChild == current_sp:
                        if len(current_sp.event_list) < 1:
                            print(3)
                            co = Tree.Tree()
                            co.leftChild = current_sp.deepcopy_single()
                            co.rightChild = new_topo.deepcopy_single()

                            self.label_lost_child(co.rightChild)

                            co.leftChild.parent = co
                            co.rightChild.parent = co

                            current_sp.parent.leftChild = co
                            current_sp.parent.isLeaf = None
                            
                            current_sp.parent.rightChild = current_sp.parent.rightChild
                            current_sp.parent.isLeaf = None
                            current_sp.parent.rightChild.parent = current_sp.parent
                            current_sp.parent.leftChild.parent = current_sp.parent

                            current_sp.parent.children = [current_sp.parent.leftChild, current_sp.parent.rightChild]
                        else:
                            print(4)
                            co = Tree.Tree()
                            co.leftChild = current_sp.parent.rightChild.deepcopy_single()
                            co.rightChild = new_topo.deepcopy_single()

                            self.label_lost_child(co.rightChild)

                            co.leftChild.parent = co
                            co.rightChild.parent = co

                            current_sp.parent.rightChild = co
                            current_sp.parent.isLeaf = None

                            current_sp.parent.leftChild = current_sp.parent.leftChild
                            current_sp.parent.isLeaf = None
                            current_sp.parent.rightChild.parent = current_sp.parent
                            current_sp.parent.leftChild.parent = current_sp.parent

                            current_sp.parent.children = [current_sp.parent.leftChild, current_sp.parent.rightChild]

                    else:
                        if len(current_sp.event_list) < 1:
                            co = Tree.Tree()
                            co.leftChild = current_sp.deepcopy_single()
                            co.rightChild = new_topo.deepcopy_single()

                            self.label_lost_child(co.rightChild)

                            co.leftChild.parent = co
                            co.rightChild.parent = co

                            current_sp.parent.rightChild = co
                            current_sp.parent.leftChild = current_sp.parent.leftChild
                            current_sp.parent.isLeaf = None
                            current_sp.parent.rightChild.parent = current_sp.parent
                            current_sp.parent.leftChild.parent = current_sp.parent
                        else:
                            co = Tree.Tree()
                            co.leftChild = current_sp.parent.leftChild.deepcopy_single()
                            co.rightChild = new_topo.deepcopy_single()

                            self.label_lost_child(co.rightChild)

                            co.leftChild.parent = co
                            co.rightChild.parent = co

                            current_sp.parent.leftChild = co
                            current_sp.parent.isLeaf = None

                            current_sp.parent.rightChild = current_sp.parent.rightChild
                            current_sp.parent.isLeaf = None
                            current_sp.parent.rightChild.parent = current_sp.parent
                            current_sp.parent.leftChild.parent = current_sp.parent

                        current_sp.parent.children = [current_sp.parent.leftChild, current_sp.parent.rightChild]

                current_sp.taxa = ''
                continue

            stack.append((current_re, current_sp.leftChild))
            stack.append((current_re, current_sp.rightChild))


    def copy_event_(self, re, sp, new_topo):
        stack = [sp]
        while stack:
            current_sp = stack.pop()
            if current_sp:
                if current_sp.id == re.id:
                    current_sp.event_list.append(re.event_list)
                else:
                    stack.append(current_sp.leftChild)
                    stack.append(current_sp.rightChild)




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
                        print('111')
                        sp.evolve='Loss'
                        return sp
                    else:
                        sp.evolve='Speciation'
                        return sp
                elif (sp.leftChild.isLeaf and sp.rightChild.isLeaf) :
                        print('2')
                        if sp.evolve==None:
                            sp.evolve= 'Speciation'
                        if len(left_.intersection(tr_))==0:
                            sp.leftChild.evolve='Loss'
                            sp.leftChild =sp.leftChild
                            self.copy_loss(tr,self.gene_tree,sp.leftChild)
                            sp.rightChild= self.reconcILS(tr,sp.rightChild,sp_copy,sp_,visited)
                            return sp
                        
                        if len(right_.intersection(tr_))==0:
                            sp.rightChild.evolve='Loss'
                            sp.rightChild =sp.rightChild
                            self.copy_loss(tr,self.gene_tree,sp.rightChild)
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
                        print(113232)
                        sp.evolve= 'Speciation'


                        if len(left_.intersection(tr_))==0:
                            sp.leftChild.evolve='Loss'
                            sp.leftChild =sp.leftChild
                            self.copy_loss(tr,self.gene_tree,sp.leftChild)
                            sp.rightChild= self.reconcILS(tr,sp.rightChild,sp_copy,sp_,visited)
                            return sp

                        print('-----------------------------------------------')
                        
                        if len(right_.intersection(tr_))==0:
                            sp.rightChild.evolve='Loss'
                            sp.rightChild= sp.rightChild
                            self.copy_loss(tr,self.gene_tree,sp.rightChild)
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

                                self.copy_event_(i[0],self.gene_tree,new_topo)

                        
                            


                            

                            

                            #visited.append(new_topo.to_newick())
                            
                            return self.reconcILS(new_topo,sp,sp_copy,sp_,visited)

                            
                        


                        

                    if  NNI_cost>=duplication_cost or  sp.isLeaf or new_multiple>=Initial_multiple_mapping:
                            print('Dups')
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
        gc.enable()
        
        while True:
            if sp:
                print('species',sp.to_newick())
                print('gene',tr.to_newick())
                tr_copy_1,tr_copy_2 = tr.deepcopy_double()
                

                sp_1 = sp.deepcopy_single()
                print('done')


                Initial_multiple_mapping = len(sp.refTo)

                print('Multiple_mapping', Initial_multiple_mapping)

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
                        print(1)
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
                                print('loss')
                                sp.leftChild.evolve='Loss'
                                self.copy_loss(tr,self.gene_tree,sp.leftChild)
                                eve.append('Loss')
                                stack.append((tr,sp.rightChild))

                            
                            elif len(right_.intersection(tr_))==0:
                                sp.rightChild.evolve='Loss'
                                eve.append('Loss')
                                print('loss1')
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
                                print('loss2')
                                self.copy_loss(tr,self.gene_tree,sp.leftChild)
                                sp.leftChild.evolve='Loss'
                                eve.append('Loss')
                            
                                stack.append((tr,sp.rightChild))
                                
                            elif len(right_.intersection(tr_))==0:
                                sp.rightChild.evolve='Loss'
                                eve.append('Loss')
                                print('loss3')
                                self.copy_loss(tr,self.gene_tree,sp.rightChild)
                                print('done')
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

                            new_topo,cost,bi_cos,child_=ILS.ILS().ILS(tr_copy_1,sp_1,sp_1,Initial_multiple_mapping,[])
                            
                            new_topo.reset()

                            recon_1 = sp_1.deepcopy_single()
                            recon_1.reset()
                            recon_1.label_internal()


                            new_gene_tree =tr_copy_2.deepcopy_single()
                            new_gene_tree.reset()
                            #new_gene_tree.leftChild = sp.parse(new_gene_tree.leftChild.to_newick())
                            #new_gene_tree.rightChild= sp.parse(new_gene_tree.rightChild.to_newick())recon_1


                            

                            recon_left,recon_right = sp.deepcopy_double()
                            
                            

                            #recon_right = sp.parse(recon_right.to_newick())
                            #recon_left= sp.parse(recon_left.to_newick())



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

                                eve+=['NNI' for i in range(Initial_multiple_mapping- cost)]
                                stack.append((new_topo,sp))

                                for i in list(child_):
                                    if type(i[1].numbered_taxa)==set:
                                        li= '-'.join(list(i[1].numbered_taxa))
                                    else:
                                        li =i[1].numbered_taxa
                                    li.replace(';','')


                            
                                    #i[0].event_list+=[li,['I' for i in range(Initial_multiple_mapping- cost)]]
                                    i[0].event_list+=[li,['I']]

                                    self.copy_event_(i[0],self.gene_tree,new_topo)                                
                            


                            

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
                                    

                                       
                                            
                                        sp.isLeaf=None
                                    
                                




                if not stack:
                    print(1)
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
    sp_string='((((A,(B,C)),D),(E,F)), (((G,(H,I)),J),(K,L)));'
    gene_tree='(((((D_4:4.05477584890832,((B_2:2.0662222193121957,(((D_9:0.8446130543699939,((B_14:0.7278290091270381,(B_9:0.4892480210277963,C_9:0.4892480210277963):0.23858098809924178):0.11678404524295583):0.0):0.557368639220968,(D_2:1.2229636089476192,A_2:1.2229636089476192):0.17901808464334268):0.664240525721234):0.0):1.9203250305987538,((E_5:3.5383328557702027,((((((((C_8:0.3951207596191984,A_8:0.39512075961919835):0.02717871181569187,((A_19:0.10113243408502237,A_17:0.10113243408502237):0.3211670373498679):0.0):1.5973409867760993):0.0,((A_18:0.1756897117023724,A_6:0.1756897117023724):1.8439507465086171):0.0):0.30500679881925974):0.0,((F_25:2.2469912800386123,E_6:2.2469912800386123):0.07765597699163695):0.0):0.6208325366119989,(F_3:2.3574175620791467,(C_3:1.318801040959761,((B_16:0.5395380854454159,(A_3:0.46158438219065884,B_3:0.4615843821906589):0.07795370325475703):0.779262955514345):0.0):1.0386165211193858):0.5880622315631014):0.5928530621279546):0.0):0.4482143941407466):0.0):0.06822859899737033):0.8950077861022718):0.0):0.0,((((B_1:0.8106737552428194,(C_1:0.7175893618222235,A_1:0.7175893618222235):0.09308439342059581):0.7909439729441123,(((C_10:0.9819174063226181,A_10:0.9819174063226181):0.3598040573558172,(((((((B_15:0.5493976865631844,C_15:0.5493976865631844):0.19858163396131828,C_13:0.7479793205245027):0.2981243783180789):0.0,((D_24:0.3677501268999852,D_1:0.3677501268999852):0.6783535719425964):0.0):0.0937700986577783):0.0,(((A_12:0.31943284517819664,C_12:0.3194328451781967):0.48899102581042064,((D_22:0.5047549466882123,((D_23:0.41961167725175946,D_21:0.41961167725175946):0.08514326943645284):0.0):0.30366892430040504):0.0):0.33144992651174254):0.0):0.20184766617807548):0.0):0.25989626450849634):0.0):1.0421537390732554,((A_11:1.21280374518546,(B_7:0.25694946656029377,((B_20:0.21352090230272383,C_7:0.21352090230272383):0.04342856425756994):0.0):0.9558542786251663):1.4309677220747272):0.0):2.3060121677504046):0.0):0.0);'
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
    re_w = readWrite.readWrite()

    if len(gene_tree)<50:

        reconcILS.reconcILS(tr,sp,sp_copy,sp,[])
        
    
        li =re_w.sp_event(sp,[])
    else:
        print('Using Iterative Function')
        li= reconcILS.iterative_reconcILS(tr,sp,sp_copy,sp,[])
        
        
    print(red.to_newick(reconcILS.gene_tree))
    print(li)

    #print('######################33')
    

    re_w = readWrite.readWrite()
 
    #li =re_w.sp_event(sp,[])
    


    sp.reset()

    tr.order_gene(sp)
    

    #print(Counter(li))
    dic={'Process':[],'Replicate':[],'Gene_tree':[],'Species_Tree':[],'Duplication':[],'NNI':[],'Loss':[]}
    
    dic['Gene_tree']+=[red.to_newick(reconcILS.gene_tree)]
    dic['Species_Tree']+=[sp_string]
        
    #print(li)
    dic= re_w.Create_pd('reconcILS',0,li,dic)


    
    df = pd.DataFrame(dic)
  
    
    
    df.to_csv(parser.output, index=False)

    dic_log ={'Gene_Tree':[tr.to_newick()],'Species_Tree':[sp_string],'Duplication_cost':[str(reconcILS.D_cost)],'NNI_cost':[str(reconcILS.I_cost)],'Loss_cost':[str(reconcILS.L_cost)]}
    df_log=pd.DataFrame(dic_log)
    df_log.to_csv(parser.output[:-4]+'_log.csv',index=False)
if __name__ == "__main__":
    
    main()
