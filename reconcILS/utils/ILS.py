import utils.Tree as Tree
import copy
import utils.Tally as Tally

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
        


    def pick_first_edge(self,child,gene_tree,tr,visited):

        if len(child)==0:
                return 
        else:
                pool={}
                tre_pool={}
                orientation={}
                
                for k in range(len(child)):
                        
                        ch1=child[k]
                        
                        gene_tree.reset()

                        ch = copy.deepcopy(ch1[0])
                        ch.reset()

                        list_tree= ch.NNI(gene_tree,ch1[2])
                        #list_tree = ch.NNI1(gene_tree,ch1)

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

                                tr.reset()
                                li[1].order_gene(tr)

                                tr.label_internal()
                                li[1].label_internal()
                                li[1].map_gene(tr)

                                number_map= len(tr.refTo)
                                tr.reset()

                                #print('##############################')
                                #print(li[1].to_newick())
                                #print(tr.to_newick())
                                #print('top_0',li[0].to_newick())
                                #print('topo_1',li[1].to_newick())
                                #print('sp',tr.to_newick())
                                #print('left_loss_cost',loss_left)
                                
                                #print('right_loss_cost',loss_right)
                                
                                #print('left_bi_cost',bi_score_left)
                                
                                #print('right_bi_cost',bi_score_right)  

                                number_map=1
                                #print((loss_score+bi_score_left+bi_score_right)*number_map)              
                                if k not in pool.keys():
                                    pool[k]= (loss_score+bi_score_left+bi_score_right)*number_map
                                    tre_pool[k]=li[1]
                                    orientation[k]=li[2]

                                else:
                                    if pool[k]<(loss_score+bi_score_left+bi_score_right)*number_map:
                                        continue
                                    else:
                                        pool[k] =(loss_score+bi_score_left+bi_score_right)*number_map
                                        tre_pool[k]=li[1]
                                        orientation[k]=li[2]
                                    
                


                min_key = min(pool, key=pool.get)

            
                M=True

                '''
                while M:

                    if tre_pool[min_key].to_newick()  in visited:
                        del pool[min_key]
                        if len(pool)==0:
                             return -1,-1,-1
                        else:
                            min_key = min(pool, key=pool.get)
                    else:
                '''
                return child[min_key],tre_pool[min_key],pool[min_key],orientation[min_key]

    def find_parent_child(self,root,child):
        if len(root.refTo)>1:
                for tre in root.refTo:
                    for tree1 in root.refTo:
                        
                        if (tree1 in tre.children):
                            if tree1==tre.leftChild:
                                ####print('match_left')
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

    def ILS(self,gene_tree,tr,sp_copy,cost,visited):
        list_tree=[]
        child=[]

        child= self.parent_child(tr,child)
        
        
        ##print(child)

        if len(child)==0 or cost<=1:
            return gene_tree,cost,-1,visited
        else:
            
            if len(child)>1:
                chil, trei , cos,orientation =self.pick_first_edge(child,gene_tree,tr,visited)

                if cos==0:
                        trei.label_internal()
                        Tally.Tally().tally_NNI(tr,trei,chil[2])
                        #return  self.ILS(trei,tr,sp_copy,cost-1,visited),cos
                        #print(chil)
                        if chil[2]=='Left':
                            if orientation=='left':
                                chii=[chil[0].leftChild,chil[0].leftChild.leftChild]
                            else:
                                chii=[chil[0].leftChild,chil[0].leftChild.rightChild]
                        else:
                            if orientation=='left':
                                    chii=[chil[0].rightChild,chil[0].rightChild.leftChild]
                            else:
                                chii=[chil[0].rightChild,chil[0].rightChild.rightChild]
                        visited.append(chii)
                        return  trei,cost-1,cos,visited
                elif cos==-1:
                        return gene_tree,0,-1,visited
                else:
                    child=[chil]


            
                 
            
            ###print(child)

            new_topo=copy.deepcopy(gene_tree)
            geneTree =copy.deepcopy(new_topo)
            geneTree.reset()
            
            
            

            

            
            for ch1 in child:
                    if ch1 in tr.visited:
                        ###print('visited')
                        continue 
                    ch = copy.deepcopy(ch1[0])
                    ch.reset()
                    
                    list_tree= ch1[0].NNI(geneTree,ch1[2])
                    #list_tree = ch.NNI1(gene_tree,ch1)
                    best_cost=cost
                    imporvement=False

                    new_topo=copy.deepcopy(geneTree)
                    
                    for i in list_tree:
                        i[1].reset()
                        i[0].reset()
                        cop= copy.deepcopy(sp_copy)
                        cop.reset()
                            
                        #print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
                        #print('top_0',i[0].to_newick())
                        #print('topo_1',i[1].to_newick())
                        #print('sp',cop.to_newick())

                        
                        #new_cost =cop.optimize_cost(i[0],i[1])
                        i[1].order_gene(cop)
                        i[1].label_internal()
                        cop.label_internal()
                        i[1].map_gene(cop)
                        new_cost = len(cop.refTo)

                        if best_cost>new_cost and cost>0:
                            if best_cost!=new_cost:
                                imporvement=True
                            else:
                                 imporvement=False
                            best_cost=new_cost


                            if tr.isRoot:
                                
                                new_topo=copy.deepcopy(i[1])
                            else:
                                new_topo=copy.deepcopy(i[1])

                            if child[0][2]=='Left':
                                if i[2]=='left':
                                    chii=[child[0][0].leftChild,child[0][0].leftChild.leftChild]
                                else:
                                    chii=[child[0][0].leftChild,child[0][0].leftChild.rightChild]
                            else:
                                if i[2]=='left':
                                    chii=[child[0][0].rightChild,child[0][0].rightChild.leftChild]
                                else:
                                    chii=[child[0][0].rightChild,child[0][0].rightChild.rightChild]
                            visited.append(chii)
                            #visited.append(i[0].to_newick())
                            cost=cost-1

                            Tally.Tally().tally_NNI(tr,new_topo,ch1[2])  

                           


                    #cost=cost-1
                    tr.visited.append([ch1[0],ch1[1]])
                    
                    if cost<=0 or imporvement==False:
                            return new_topo,cost,-1,visited
                    else:
                            new_sp = copy.deepcopy(sp_copy)
                            new_sp.reset()
                            new_topo.reset()
                            new_topo.order_gene(new_sp)
                            
                            new_topo.label_internal()
                        
                            new_sp.label_internal()
                            
                            new_topo.map_gene(new_sp)

                            

                            return self.ILS(new_topo,new_sp,sp_copy,cost,visited)
                    
                    
        return new_topo,cost,-1,visited
