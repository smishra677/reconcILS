from .utils_reconcILS import Tree
import copy
from .utils_reconcILS import Tally

class ILS:
    def find_bipartitions(self,bipartitions, subtree):
        if not subtree:
            return []

        stack = [subtree]

        while stack:
            node = stack.pop()

            if not node.isLeaf:
                bipartitions.append(node.taxa_list)
                
                if node.rightChild:
                    stack.append(node.rightChild)
                if node.leftChild:
                    stack.append(node.leftChild)

        return bipartitions

    def find_biparition_cost(self,sp,tr):

        
        merged = sp+tr



        difference= 0

        diff_parition=[]
        for i in merged:
            if i in tr and i not in sp:
                diff_parition.append(i)
                difference= difference+1
    
        
        return difference
        

    def find_number_map(self,sp):
        stack=[sp]
        total_map=0
        while stack:
            current_node=stack.pop()

            if current_node: 
                map_=len(current_node.refTo)
                if map_>1:
                    total_map+=map_

                stack.append(current_node.leftChild)

                stack.append(current_node.rightChild)
        
        return total_map


    def find_lowest(self,gene_tree,id_list,child_dic):
        

        level_dic={}

        address_list= {child_dic[id][1]:id for id in id_list}

        #print(id_list)
        #print(address_list)
        for add in address_list.keys():
            stack=[(add,0)]
            while stack:
                current_node,level=stack.pop()
                
                if current_node:
                    if current_node.isLeaf:
                        
                        level_dic[address_list[add]]=level

                    stack.append((current_node.leftChild,level+1))
                    stack.append((current_node.rightChild,level+1))
                
        id_max_level= min(level_dic, key=level_dic.get)
        #print(level_dic)
        #print([i.to_newick() for i in address_list.keys()])

        return id_max_level

    
    def pick_first_edge(self,child,gene_tree,tr,visited):

            if len(child)==0:
                return 
            else:
                pool={}
                tre_pool={}
                orientation={}
                super_list={}
                LCA_dic={}
                for k in range(len(child)):
                        
                        ch1=child[k]
                        
                        gene_tree.reset()

                        ch = ch1[0].deepcopy()
                        ch.reset()


                        list_tree= ch.NNI(gene_tree,ch1[2])

                        super_list[k]=list_tree
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
                                number_map_left=len(tr.refTo)
                                



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

                                #number_map= len(tr.refTo)
                                number_map_right=len(tr.refTo)
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
                                #firstTree = ete3.Tree(tr.to_newick())
                                #secondTree = ete3.Tree(li[1].to_newick())
                                #rf, _, _, _, _, _, _ = firstTree.robinson_foulds(secondTree)

                                #number_map=1
                                #print((loss_score+bi_score_left+bi_score_right)*number_map)              
                                if number_map_right==0:
                                    number_map_right=1
                                if number_map_left==0:
                                    number_map_left=1

                                if k in LCA_dic:
                                    LCA_dic[k] = min(LCA_dic[k],number_map_left+number_map_right)
                                else:
                                    LCA_dic[k] = number_map_left+number_map_right

                                combined_score =(loss_left+bi_score_left)+(loss_right+bi_score_right)
                                
                                

                                if k not in pool:
                                    pool[k] = combined_score
                                    tre_pool[k] = li[1]
                                    orientation[k] = li[2]
                                else:
                                    if pool[k] > combined_score:
                                        pool[k] = combined_score
                                        tre_pool[k] = li[1]
                                        orientation[k] = li[2]



                min_key = min(pool, key=pool.get)

                
                
                min_combined_score=min(pool.values())
                min_value_keys_count = sum(1 for value in pool.values() if value == min_combined_score)


                
                if min_value_keys_count>1:
                    keys_with_min_value = [key for key, value in pool.items() if value == min_combined_score]
                    extracted_dic = {key: LCA_dic[key] for key in keys_with_min_value}
                    min_key = min(extracted_dic, key=extracted_dic.get)
                    min_combined_score_1=min(extracted_dic.values())
                    min_value_keys_count_1 = sum(1 for value in extracted_dic.values() if value == min_combined_score_1)
                    keys_with_min_value_1 = [key for key, value in extracted_dic.items() if value == min_combined_score_1]
                    if min_value_keys_count_1>1:
                        min_key=self.find_lowest(gene_tree,keys_with_min_value_1,
                        child)

                
                        
                

                

                        

                #print("==============================<><()))()()()()()",pool,LCA_dic)

                return child[min_key],tre_pool[min_key],pool[min_key],orientation[min_key],super_list[min_key]
        


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

    def ILS(self, gene_tree, tr, sp_copy, cost,best_cost, visited):
        #print(1,cost)
        #Initial_best_cost=best_cost
        child = self.parent_child(tr, [])
        same_round=False
        if len(child) == 0 or cost <= 0:
            return gene_tree, cost, -1, visited

        geneTree = gene_tree.deepcopy()
        geneTree.reset()
        
        if len(child) == 1:
            list_tree = child[0][0].NNI(geneTree, child[0][2])
        else:
            chil, trei, cos, orientation, list_tree = self.pick_first_edge(child, gene_tree, tr, visited)
            #print(cos,trei.to_newick())
            if cos == 0:
                trei.label_internal()
                chii = self.get_child_info(chil, orientation)
                visited.append(chii)
                return trei, cost - 1, cos, visited
            elif cos == -1:
                return gene_tree, 0, -1, visited
            else:
                child = [chil]

        ch = child[0][0].deepcopy()
        ch.reset()

        #best_cost = cost
        improvement = False

        new_topo = geneTree.deepcopy()
        
        for i in list_tree:
            i[1].reset()
            i[0].reset()
            cop = sp_copy.deepcopy()
            cop.reset()

            i[1].order_gene(cop)
            i[1].label_internal()
            cop.label_internal()
            i[1].map_gene(cop)
            new_cost = len(cop.refTo)

            if best_cost >= new_cost and cost >= 0:
                improvement = best_cost >= new_cost
                best_cost = new_cost
                new_topo = i[1].deepcopy()

                chii = self.get_child_info(child[0], i[2])
                
                

 
                if same_round:
                    visited=visited[:-1]
                else:
                    same_round=True
                    cost -=1
                visited.append(chii)
                
            
            
        if cost <= 0 or not improvement:
            return new_topo, cost, -1, visited

        new_sp = sp_copy.deepcopy()
        new_sp.reset()
        new_topo.reset()
        new_topo.order_gene(new_sp)
        new_topo.label_internal()
        new_sp.label_internal()
        new_topo.map_gene(new_sp)

        return self.ILS(new_topo, new_sp, sp_copy, cost,best_cost, visited)

    def get_child_info(self, chil, orientation):
        if chil[2] == 'Left':
            if orientation == 'left':
                return [chil[0].leftChild, chil[0].leftChild.leftChild]
            else:
                return [chil[0].leftChild, chil[0].leftChild.rightChild]
        else:
            if orientation == 'left':
                return [chil[0].rightChild, chil[0].rightChild.leftChild]
            else:
                return [chil[0].rightChild, chil[0].rightChild.rightChild]
