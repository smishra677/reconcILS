import utils.Tree as Tree
import re
import os
import copy
import utils.Tally as Tally
import argparse

class ILS:
    def find_bipartitions(self,bi,sp):

        if sp:
            if sp.isLeaf:
                return bi
            else:
                bi.append(sp.taxa_list)
            bi1 =self.find_bipartitions(bi,sp.leftChild)
            bi2 =self.find_bipartitions(bi1,sp.rightChild)
        return bi2

    def find_biparition_cost(self,sp,tr):

        
        merged = sp+tr



        difference= 0

        diff_parition=[]
        for i in merged:
            if i in tr and i not in sp:
                diff_parition.append(i)
                difference= difference+1
    
        
        return difference
        


    def pick_first_edge(self,child,gene_tree,tr):

        if len(child)==0:
                return 
        else:
                pool={}
                tre_pool={}
                
                for k in range(len(child)):
                        
                        ch1=child[k]
                        
                        gene_tree.reset()

                        ch = copy.deepcopy(ch1[0])
                        ch.reset()

                        list_tree= ch.NNI(gene_tree,ch1[2])


                        for li in list_tree:
                            
                               
                                
                                tr.reset()
                                li[0].reset()
                                li[0].order_gene(tr)

                                li[0].label_internal()
                                tr.label_internal()
                                bi_sp= self.find_bipartitions([],tr)
                                bi_li_0_l= self.find_bipartitions([],li[0].leftChild)
                                bi_li_0_r= self.find_bipartitions([],li[0].rightChild)
                                bi_score_left= self.find_biparition_cost(bi_sp,bi_li_0_l)
                                bi_score_right= self.find_biparition_cost(bi_sp,bi_li_0_r)

                                tr.reset()
                                tr.cost=0
                                li[0].reset()

                                li[0].leftChild.order_gene(tr)
                                li[0].leftChild.label_internal()
                                tr.label_internal()
    
                                li[0].leftChild.map_gene(tr)
                                tr.find_loss_sp(tr)
                                loss_left = tr.cost


                                
                                tr.reset()
                                tr.cost=0
                                li[0].rightChild.order_gene(tr)
                                li[0].rightChild.label_internal()
                                tr.label_internal()

                                
                                li[0].rightChild.map_gene(tr)
                                tr.find_loss_sp(tr)
                                loss_right = tr.cost
                                loss_score=loss_left+loss_right

                        
                                if k not in pool.keys():
                                    pool[k]= loss_score+bi_score_left+bi_score_right
                                    tre_pool[k]=li[1]

                                else:
                                    if pool[k]<(loss_score+bi_score_left+bi_score_right):
                                        continue
                                    else:
                                        pool[k] =loss_score+bi_score_left+bi_score_right
                                        tre_pool[k]=li[1]
                                    

                return child[min(pool, key=pool.get)],tre_pool[min(pool, key=pool.get)],pool[k]

    def find_parent_child(self,root,child):

        if len(root.refTo)>1:
                #root.refTo.reverse()

                for tre in root.refTo:
                    for tree1 in root.refTo:
                        
                        if (tree1 in tre.children):
                            if tree1==tre.leftChild:
                                #print('match_left')
                                child.append([tre,tree1,'Left'])
                            else:
                                child.append([tre,tree1,'Right'])

        return child
    

    def parent_child(self,root,child):
        if root:
            if root.isLeaf:
                return []
            else:
                child= self.find_parent_child(root,child)
        return child

    def ILS(self,gene_tree,tr,sp_copy,cost):
        list_tree=[]
        child=[]

        child= self.parent_child(tr,child)

    
        if len(child)==0 or cost==0:
            return gene_tree,cost
        else:
            
            if len(child)>1:
                chil, trei , cos =self.pick_first_edge(child,gene_tree,tr)
                
                if cos==0:
                        trei.label_internal()
                        Tally.Tally().tally_NNI(tr,trei,chil[2])
                        return  self.ILS(trei,tr,sp_copy,cost-1)
                else:
                    child=[chil]



           
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
                            


                        
                        #new_cost =cop.optimize_cost(i[0],i[1])
                        i[1].order_gene(cop)
                        i[1].label_internal()
                        cop.label_internal()
                        i[1].map_gene(cop)
                        new_cost = len(cop.refTo)

                        if best_cost>new_cost and cost>0:
                            best_cost=new_cost


                            if tr.isRoot:
                                
                                new_topo=copy.deepcopy(i[1])
                            else:
                                new_topo=copy.deepcopy(i[1])
                            
                            print(88888888888888888888888888)
                            Tally.Tally().tally_NNI(tr,new_topo,ch1[2])  
                    
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

                            

                            return self.ILS(new_topo,new_sp,sp_copy,cost)
                    
                    
        return new_topo,cost
