import utils.Tree as Tree
import copy
import utils.Tally as Tally

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
        


    def pick_first_edge(self,child,gene_tree,tr,visited):

        if len(child)==0:
                return 
        else:
                pool={}
                tre_pool={}
                orientation={}
                super_list={}
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

                                print('##############################')
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

                                combined_score = (loss_score + bi_score_left + bi_score_right) * number_map

                                if k not in pool:
                                    pool[k] = combined_score
                                    tre_pool[k] = li[1]
                                    orientation[k] = li[2]
                                else:
                                    if pool[k] >= combined_score:
                                        pool[k] = combined_score
                                        tre_pool[k] = li[1]
                                        orientation[k] = li[2]



                min_key = min(pool, key=pool.get)


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

    def ILS(self, gene_tree, tr, sp_copy, cost, visited):
        child = self.parent_child(tr, [])

        if len(child) == 0 or cost <= 1:
            return gene_tree, cost, -1, visited

        geneTree = gene_tree.deepcopy()
        geneTree.reset()
        
        if len(child) == 1:
            list_tree = child[0][0].NNI(geneTree, child[0][2])
        else:
            chil, trei, cos, orientation, list_tree = self.pick_first_edge(child, gene_tree, tr, visited)

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

        best_cost = cost
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

            if best_cost > new_cost and cost > 0:
                improvement = best_cost != new_cost
                best_cost = new_cost
                new_topo = i[1].deepcopy()

                chii = self.get_child_info(child[0], i[2])
                visited.append(chii)

                cost -= 1

        if cost <= 0 or not improvement:
            return new_topo, cost, -1, visited

        new_sp = sp_copy.deepcopy()
        new_sp.reset()
        new_topo.reset()
        new_topo.order_gene(new_sp)
        new_topo.label_internal()
        new_sp.label_internal()
        new_topo.map_gene(new_sp)

        return self.ILS(new_topo, new_sp, sp_copy, cost, visited)

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
