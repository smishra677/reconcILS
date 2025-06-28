from itertools import count
import copy
import matplotlib.pyplot as plt
from . import idmaker_
from . import readWrite
import pickle


class Tree:
    def __init__(self):
        self.id= idmaker_.idmaker2().id
        self.taxa=None
        self.numbered_taxa=None
        self.event=None
        self.evolve=None
        self.tag = None
        self.isRoot= None
        self.isLeaf =None
        self.leftChild=None
        self.rightChild=None
        self.children=[]
        self.refTo=[]
        self.parent=None
        self.split_list=None
        self.to_tag=None
        self.cost=0
        self.inital_ref=0
        self.event_list=[]
        #self.visited=[]
        #self.li=[]
        #self.sp_ev_list=[]
        #self.paralogy=[]
        #self.NNI_=[]
        self.taxa_list=[]

    
    def deepcopy_single(self):
        #deep_copied_instance = pickle.loads(pickle.dumps(self))
    
        return pickle.loads(pickle.dumps(self))

    
    def deepcopy(self):
        #deep_copied_instance = pickle.loads(pickle.dumps(self))

        return pickle.loads(pickle.dumps(self))

    
    def deepcopy_double(self):
        deep_copied_instance = pickle.dumps(self)
        copy_1=pickle.loads(deep_copied_instance)
        copy_2=pickle.loads(deep_copied_instance)
        del deep_copied_instance

        return copy_1,copy_2


    
    def clear_ref(self):
        stack = [self]
        while stack:
            node = stack.pop()
            if node:
                node.refTo = []
                node.tag = None
                stack.append(node.leftChild)
                stack.append(node.rightChild)

    
    def reset(self):
        stack = [self]
        while stack:
            node = stack.pop()
            if node:
                node.event = None
                node.tag = None
                node.evolve = None
                node.taxa_list = []
                node.cost=0
                #node.sp_ev_list=[]
                #node.NNI_ = []

                if node.isLeaf == None:
                    node.taxa = None
                node.refTo = []
                stack.append(node.leftChild)
     
                stack.append(node.rightChild)

    '''
    def reset_tag(self):
        stack = [self]
        while stack:
            node = stack.pop()
            if node:
                node.tag = None
                stack.append(node.leftChild)
                stack.append(node.rightChild)





    def reset_evolve(self):
        stack = [self]
        while stack:
            node = stack.pop()
            if node:
                if node.evolve == 'Loss':
                    node.evolve = None
                stack.append(node.leftChild)
                stack.append(node.rightChild)

    '''

    
    def printorder_species(self, root):
        stack = [(root, False)]
        while stack:
            node, visited = stack.pop()
            if node:
                if visited:
                    if self.isLeaf and (self.taxa == node.taxa) and self.tag is None:
                        self.tag = node
                        node.refTo.append(self)
                else:
                    stack.append((node, True))
                    stack.append((node.rightChild, False))
                    stack.append((node.leftChild, False))

    
    def search_sp(self,species_tree):
        self.printorder_species(species_tree)
        pass
    
    
    def order_gene(self, species_tree):
        stack = [self]
        
        while stack:
            current_node = stack.pop()

            if current_node:
                current_node.search_sp(species_tree)

                if current_node.rightChild:
                    stack.append(current_node.rightChild)
                if current_node.leftChild:
                    stack.append(current_node.leftChild)

    
    def label_internal(self):
        if self != None:
            if self.leftChild:
                self.leftChild.label_internal()
                if self.tag==None:
                    self.tag= [self.leftChild.tag]
                    if type(self.leftChild.taxa)==set:
                        self.taxa=set()
                        self.taxa=self.taxa.union(self.leftChild.taxa)
                    else:
                        self.taxa=set()
                        self.taxa.add(self.leftChild.taxa)


                    if type(self.leftChild.numbered_taxa)==set:
                        self.numbered_taxa=set()
                        self.numbered_taxa= self.numbered_taxa.union(self.leftChild.numbered_taxa)
                    else:
                        self.numbered_taxa=set()
                        self.numbered_taxa.add(self.leftChild.numbered_taxa)
                    
                    if self.leftChild.isLeaf:
                        self.taxa_list.append(self.leftChild.taxa)
                        self.taxa_list.sort()
                    else:
                        self.leftChild.taxa_list.sort()
                        self.taxa_list.append(''.join(self.leftChild.taxa_list))
                        self.taxa_list.sort()

            if  self.rightChild:
                self.rightChild.label_internal()
                if  type(self.tag)==list:
                    self.tag+= [self.rightChild.tag]
                    #self.tag = list(set(self.tag))
                    if type(self.rightChild.taxa)==set:
                        #print(self.taxa,self.rightChild.taxa)
                        self.taxa=self.taxa.union(self.rightChild.taxa)
                        #print(self.taxa)
                    else:
                        self.taxa.add(self.rightChild.taxa)
                    if type(self.rightChild.numbered_taxa)==set:
                        self.numbered_taxa= self.numbered_taxa.union(self.rightChild.numbered_taxa)
                    else:
                        self.numbered_taxa.add(self.rightChild.numbered_taxa)


                    if self.rightChild.isLeaf:
                        self.taxa_list.append(self.rightChild.taxa)
                        self.taxa_list.sort()

                    else:
                        self.rightChild.taxa_list.sort()
                        self.taxa_list.append(''.join(self.rightChild.taxa_list))  
                        self.taxa_list.sort()               

    

    
    def map_species(self, root):
        stack = [(root, False)]
        while stack:
            node, visited = stack.pop()
            if node:
                if visited:
                    if self.isLeaf is None and node.isLeaf is None and len(self.taxa) > 1:
                        if set(self.taxa).issubset(node.taxa) and self.event is None:
                            self.event = node
                            node.refTo.append(self)
                else:
                    stack.append((node, True))
                    stack.append((node.rightChild, False))
                    stack.append((node.leftChild, False))


                    
    
    def search_spMap(self,species_tree):
        self.map_species(species_tree)
        pass


    
    def map_gene(self, species_tree):
        stack = [(self, False)]
        while stack:
            node, visited = stack.pop()
            if node:
                if visited:
                    node.search_spMap(species_tree)
                else:
                    stack.append((node, True))
                    if node.rightChild:
                        stack.append((node.rightChild, False))
                    if node.leftChild:
                        stack.append((node.leftChild, False))


    
    
    def check_below(self):
        stack = [self]
        flag = False
        
        while stack:
            current_node = stack.pop()
            if current_node:
                if len(current_node.refTo) > 0:
                    flag = True
                    break
                
                if current_node.leftChild:
                    stack.append(current_node.leftChild)
                if current_node.rightChild:
                    stack.append(current_node.rightChild)
        
        return 1 if not flag else 0

    def find_loss_sp(self, root):
        if self:
            
            
            if self.leftChild and self.rightChild:
                if self.rightChild.isLeaf and len(self.rightChild.refTo) == 0 and self.leftChild.isLeaf and len(self.leftChild.refTo) == 0:
                    #print(2, root.cost)
                    root.cost+=1
                    return 0
            #if len(self.refTo)==0:
            #print('tag', self.to_newick())
            check= self.check_below()
            #print(check)
            if check==1:
                root.cost+=1
                return 0
            
            
            else:
                         
                if self.leftChild:
                    #print(3, root.cost)
                    value_left=self.leftChild.find_loss_sp(root)
                    root.cost += value_left
                    #print('tag', self.to_newick())
                    #print('ack', 3)
                    #print('ack', 3, root.cost)
                
                if self.rightChild:
                    #print(4, root.cost)
                    value_right=self.rightChild.find_loss_sp(root)
                    root.cost += value_right
                    #print('tag', self.to_newick())
                    #print('ack', 4, root.cost)
                return 0
                




    
    def change_id(self, tree, new):
        stack = [(tree, new)]
        while stack:
            tree_node, new_node = stack.pop()
            if tree_node and new_node:
                tree_node.id = new_node.id
                val_left = val_right = 0

                if tree_node.leftChild and new_node.leftChild:
                    val_left = len(set(tree_node.leftChild.taxa).intersection(set(new_node.leftChild.taxa)))
                if tree_node.rightChild and new_node.rightChild:
                    val_right = len(set(tree_node.rightChild.taxa).intersection(set(new_node.rightChild.taxa)))

                if val_left > val_right:
                    stack.append((tree_node.rightChild, new_node.rightChild))
                    stack.append((tree_node.leftChild, new_node.leftChild))
                else:
                    stack.append((tree_node.rightChild, new_node.leftChild))
                    stack.append((tree_node.leftChild, new_node.rightChild))

           
    
    def write_events(self):
        from collections import Counter
        if self.isLeaf:
            
            dic=dict(Counter(self.event_list))
            val= [str(k) + ':' + str(v) for k,v in dic.items()]
            if len(val)>0:
                return self.taxa+':'+"["+','.join(val)+"]"
            else:
                return self.taxa
            
        else:
            dic=dict(Counter(self.event_list))
            val= [str(k) + ':' + str(v) for k,v in dic.items()]
            val1= "["+','.join(val)+"]"
            return val1 if len(val)>0 else  ''


    # Got it From StackOverflow:
    #https://stackoverflow.com/questions/61117131/how-to-convert-a-binary-tree-to-a-newick-tree-using-python
    # https://stackoverflow.com/questions/61117131/how-to-convert-a-binary-tree-to-a-newick-tree-using-python
    
    def traverse(self, newick):
        if self.leftChild and not self.rightChild:
            newick = f"(,{self.leftChild.traverse(newick)}){self.taxa if self.isLeaf else ''}"
        elif not self.leftChild and self.rightChild:
            newick = f"({self.rightChild.traverse(newick)},){self.taxa if self.isLeaf else ''}"
        elif self.leftChild and self.rightChild:
            newick = f"({self.rightChild.traverse(newick)},{self.leftChild.traverse(newick)}){self.taxa if self.isLeaf else ''}"
        elif not self.leftChild and not self.rightChild :
            newick = f"{self.taxa if self.isLeaf else ''}"
        else:
            pass
        return newick





    # Got it From StackOverflow:
    #https://stackoverflow.com/questions/61117131/how-to-convert-a-binary-tree-to-a-newick-tree-using-python
    # https://stackoverflow.com/questions/61117131/how-to-convert-a-binary-tree-to-a-newick-tree-using-python
    
    def to_newick(self):
        newick = ""
        newick = self.traverse(newick)
        newick = f"{newick};"
        return newick


    #Got the framework of the code from StackOver FLow 
    # https://stackoverflow.com/questions/51373300/how-to-convert-newick-tree-format-to-a-tree-like-hierarchical-object
    # https://stackoverflow.com/questions/51373300/how-to-convert-newick-tree-format-to-a-tree-like-hierarchical-object
    
    def parse(self,newick):
        import re

        tokens = re.finditer(r"([^:;,()\s]*)(?:\s*:\s*([\d.]+)\s*)?([,);])|(\S)", newick+";")

        def recurse(tre,nextid = 0, parentid = -1):
            thisid = nextid
            #tre.id= nextid
            children = []

            name, length, delim, ch = next(tokens).groups(0)
            #print(name)
            tre.taxa= name
            if name!=0:
                tre.taxa= ''.join(name.split('_')[:2])
                tre.numbered_taxa=name
                tre.isLeaf= True
            if ch == "(":
                while ch in "(,":
                    new_tre= Tree()
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

        tre= Tree()
        val =recurse(tre)

        return val[-1]
    

    
    def find_address(self,node):
        if self:
            if self.id==node.id:
                return self
            if self.leftChild:
                return self.leftChild.find_address(node)
            if self.rightChild:
                return self.rightChild.find_address(node)

  
        
    
    def locate_copy(self, copy_tree, tree):
        stack = [tree]
        while stack:
            node = stack.pop()
            if node:
                if self.id == node.id:
                    if node.parent.leftChild == node:
                        node.parent.children.remove(node.parent.leftChild)
                        node.parent.leftChild = copy_tree
                        node.leftChild = copy_tree.leftChild
                        node.rightChild = copy_tree.rightChild
                        copy_tree.parent = node.parent
                        node.parent.children.append(copy_tree)
                        return
                    else:
                        node.parent.children.remove(node.parent.rightChild)
                        node.parent.rightChild = copy_tree
                        node.leftChild = copy_tree.leftChild
                        node.rightChild = copy_tree.rightChild
                        copy_tree.parent = node.parent
                        node.parent.children.append(copy_tree)
                        return
                else:
                    stack.append(node.rightChild)
                    stack.append(node.leftChild)        


    
    def NNI(self,gene_tree,flag):
        red = readWrite.readWrite()

        geneTree_left =gene_tree.deepcopy()
        geneTree_left.reset()

        



        geneTree_right =gene_tree.deepcopy()
        geneTree_right.reset()

        copy_left = self.deepcopy()
        copy_left.reset()
        copy_right = self.deepcopy()
        copy_right.reset()
        if flag=='Left':
            copy_right_child= self.leftChild.rightChild.deepcopy()
            copy_right_child.reset()
            copy_left_child= self.leftChild.leftChild.deepcopy()
            copy_left_child.reset()


            new_tree_left = Tree()
            new_tree_left.leftChild =copy_left_child
            new_tree_left.rightChild=copy_left.rightChild
            new_tree_left.children =[new_tree_left.leftChild, new_tree_left.rightChild]
            new_tree_left.leftChild.parent=new_tree_left
            new_tree_left.rightChild.parent=new_tree_left
            new_tree_left.id =self.leftChild.id

            


            copy_left.rightChild=new_tree_left
            copy_left.leftChild=copy_right_child
            copy_left.children= [copy_left.rightChild, copy_left.leftChild]
            new_tree_left.parent=copy_left
            copy_right_child.parent=copy_left

            new_tree_right = Tree()
            new_tree_right.leftChild=copy_right_child
            new_tree_right.rightChild=copy_right.rightChild
            new_tree_right.children = [new_tree_right.rightChild, new_tree_right.leftChild]
            new_tree_right.rightChild.parent=new_tree_right
            new_tree_right.leftChild.parent=new_tree_right


            new_tree_right.id =self.leftChild.id

            copy_right.rightChild=new_tree_right
            copy_right.leftChild=copy_left_child
            copy_right.children= [copy_right.rightChild, copy_right.leftChild]
            copy_right.rightChild.parent=copy_right
            copy_right.leftChild.parent=copy_right



        
        else:
            copy_right_child= self.rightChild.rightChild.deepcopy()
           
            copy_right_child.reset()
            copy_left_child= self.rightChild.leftChild.deepcopy()
            copy_left_child.reset()
            
            new_tree_left = Tree()
            new_tree_left.rightChild =copy_left_child
            new_tree_left.leftChild=copy_left.leftChild

            new_tree_left.children =[new_tree_left.leftChild, new_tree_left.rightChild]
            new_tree_left.leftChild.parent=new_tree_left
            new_tree_left.rightChild.parent=new_tree_left

            new_tree_left.id =self.rightChild.id



            copy_left.rightChild=new_tree_left
            copy_left.leftChild=copy_right_child
            copy_left.children= [copy_left.rightChild, copy_left.leftChild]
            copy_left.rightChild.parent=copy_left
            copy_left.leftChild.parent=copy_left

            
            new_tree_right = Tree()
            new_tree_right.leftChild=copy_right_child
            new_tree_right.rightChild=copy_right.leftChild
            #print(new_tree_right.to_newick())
            new_tree_right.children= [new_tree_right.leftChild, new_tree_right.leftChild]
            new_tree_right.leftChild.parent=new_tree_right
            new_tree_right.leftChild.parent=new_tree_right

            #copy_right_child.id=new_tree_right.id
            new_tree_right.id =self.rightChild.id


            copy_right.rightChild=new_tree_right
            copy_right.leftChild=copy_left_child
            copy_right.children= [copy_right.rightChild, copy_right.leftChild]
            copy_right.rightChild.parent=copy_right
            copy_right.leftChild.parent=copy_right        



        if self.parent==None:

            l_t= red.parse(red.to_newick(copy_left))
            r_t=red.parse(red.to_newick(copy_right))


            l_t.label_internal()
            r_t.label_internal()
            copy_left.label_internal()
            copy_right.label_internal()

            
            l_t.change_id(l_t,copy_left)
            r_t.change_id(r_t,copy_right)
            
            r_t.reset()
            copy_right.reset()
            
            l_t.reset()
            copy_left.reset()
            


            return [[l_t,l_t,'left'],[r_t,r_t,'right']]
        
        geneTree_left.label_internal()
        
        geneTree_right.label_internal()

        copy_left.label_internal()
        copy_right.label_internal()


        
        self.locate_copy(copy_left,geneTree_left)
        self.locate_copy(copy_right,geneTree_right)


        
        l_t= red.parse(red.to_newick(copy_left))
        l_t_1=red.parse(red.to_newick(geneTree_left))

        


        r_t= red.parse(red.to_newick(copy_right))
        r_t_1=red.parse(red.to_newick(geneTree_right))

        


        l_t.label_internal()
        r_t.label_internal()
        r_t_1.label_internal()
        l_t_1.label_internal()

        copy_left.label_internal()
        copy_right.label_internal()
        geneTree_left.label_internal()
        geneTree_right.label_internal()

        
        
        
                
            
            
        l_t.change_id(l_t,copy_left)
        r_t.change_id(r_t,copy_right)
        l_t_1.change_id(l_t_1,geneTree_left)
        r_t_1.change_id(r_t_1,geneTree_right)
        

        r_t.reset()
        copy_right.reset()
        l_t_1.reset()
        geneTree_left.reset()
        r_t_1.reset()
        geneTree_right.reset()
        l_t.reset()
        copy_left.reset()
        
        return [[l_t,l_t_1,'left'],[r_t,r_t_1,'right']]