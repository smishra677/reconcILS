from itertools import count
import copy
class Tree:
    new_id = count()
    def __init__(self):
        self.id= next(Tree.new_id)
        self.taxa=None
        self.event=None
        self.tag = None
        self.isRoot= None
        self.isLeaf =None
        self.leftChild=None
        self.rightChild=None
        self.evolve=None
        self.children=[]
        self.refTo=[]
        self.parent=None


    def reset(self):
        if self:
            self.event=None
            self.tag = None
            self.evolve=None
            if self.isLeaf==None:
                self.taxa=None
            self.refTo=[]

            if self.leftChild:
                self.leftChild.reset()
            
            if self.rightChild:
                self.rightChild.reset()


    
    def __setNewick__(self,newick):
        self.taxa=None


    def __setChild__(self,leftChild,rightChild):
        self.leftChild=leftChild
        self.rightChild=rightChild

    
    def __setRoot__(self):
        self.isRoot=True
        self.isLeaf=False

    def __setLeaf__(self):
        self.isLeaf=True
        self.isRoot=False


    def printorder_species(self,root):
        if root:
            if self.isLeaf==True and (self.taxa== root.taxa) and self.tag == None:
                self.tag= root
                root.refTo.append(self)
            else:
                self.printorder_species(root.leftChild)
                self.printorder_species(root.rightChild)
                

    def search_sp(self,species_tree):
        self.printorder_species(species_tree)
        pass
    

    def printorder_gene(self,species_tree):
    
        if self != None:
            self.search_sp(species_tree),
            if self.leftChild:
                self.leftChild.printorder_gene(species_tree)
            if  self.rightChild:
                self.rightChild.printorder_gene(species_tree)
    



    def label_internal(self):
        if self != None:
            if self.leftChild:
                self.leftChild.label_internal()
                if self.tag==None:
                    self.tag= [self.leftChild.tag]
                    self.taxa=set(self.leftChild.taxa)

            if  self.rightChild:
                self.rightChild.label_internal()
                if  type(self.tag)==list:
                    self.tag+= [self.rightChild.tag]
                    #self.tag = list(set(self.tag))
                    self.taxa=self.taxa.union(self.rightChild.taxa)

    


    def map_species(self,root):
        if root and  self.isLeaf==None and root.isLeaf==None:
            self.map_species(root.leftChild)
            self.map_species(root.rightChild)
            #print(self.taxa)
            #print(root.taxa)
            #print(len(self.taxa.difference(root.taxa))==0)
            #print(self.event)
            if len(self.taxa.difference(root.taxa))==0 and self.event == None:
                self.event= root
                root.refTo.append(self)



                    

    def search_spMap(self,species_tree):
        self.map_species(species_tree)
        pass

    def printorder(self):
    
        if self:
                print(self.taxa),
                print(self.event),
                print('###################'),
                if self.leftChild:
                    self.leftChild.printorder()
                if self.rightChild:
                    self.rightChild.printorder()


    def map_gene(self,species_tree):
        if self != None:
            self.search_spMap(species_tree),
            if self.leftChild:
                self.leftChild.map_gene(species_tree)
            if  self.rightChild:
                self.rightChild.map_gene(species_tree)
    

        pass


    def tag_species(self):
        if self != None:
            if len(self.refTo)>1:
                self.event='D'

            if self.leftChild:
                self.leftChild.tag_species()

            if self.rightChild:
                self.rightChild.tag_species()

           
            
    

        pass


    def tag_gene(self):
        if self != None:
            if len(self.leftChild.taxa & self.rightChild.taxa  )>1:
                self.event='D'
            else:
                self.event='S'

            if self.leftChild:
                self.leftChild.tag_gene()

            if self.rightChild:
                self.rightChild.tag_gene()

           
            
    

        pass


    

    



    

    def sp_tag(self):
    
        if self:
            if self.leftChild:
                self.leftChild.sp_tag()
            print(self.taxa),
            print(self.refTo),
            #print(root.isLeaf),
            if self.rightChild:
                self.rightChild.sp_tag()



    def tag_loss(self):
        if self != None:
            if self.isLeaf:
                if self.event==None and self.evolve==None:
                    self.evolve='L'
            else:
                if self.leftChild:
                    self.leftChild.tag_loss()
                if self.rightChild:
                    self.rightChild.tag_loss()
        pass

    def find_loss(self,recon_tree):
        recon_tree.map_gene(self)
        recon_tree.tag_loss()
        pass

    def find_cost(self,node):
        if self:
            if self.taxa==node.taxa:
                return len(self.refTo)
            else:
                if self.leftChild:
                    self.leftChild.find_cost(node)
                if self.rightChild:
                    self.rightChild.find_cost(node)
        pass

    def optimize_cost(self,node,tree2):
        tree2.printorder_gene(self)
        
        tree2.label_internal()
    
        self.label_internal()
        
        tree2.map_gene(self)
        
        


        print(self.find_cost(node))
        pass

    
    def NNI(self,flag):
        copy_left = copy.deepcopy(self)
        copy_left.reset()
        copy_right = copy.deepcopy(self)
        copy_right.reset()
        if flag=='left':
            copy_right_child= copy.deepcopy(self.leftChild.rightChild)
            copy_left_child= copy.deepcopy(self.leftChild.leftChild)
            new_tree_left = Tree()
            new_tree_left.leftChild =copy_left_child
            new_tree_left.rightChild=copy_left.rightChild

            copy_left.rightChild=new_tree_left
            copy_left.leftChild=copy_right_child


            new_tree_right = Tree()
            new_tree_right.leftChild=copy_right_child
            new_tree_right.rightChild=copy_right.rightChild

            copy_right.rightChild=new_tree_right
            copy_right.leftChild=copy_left_child


        
        return [copy_left,copy_right]

