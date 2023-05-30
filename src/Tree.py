class Tree:
    def __init__(self):
        self.id=None
        self.taxa=None
        self.event=None
        self.tag = None
        self.isRoot= None
        self.isLeaf =None
        self.leftChild=None
        self.rightChild=None
        self.children=[]

    
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
                print(root.taxa)
                print(self.taxa)
                self.tag= root
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
                    self.taxa=self.taxa.union(self.rightChild.taxa)

    



    def map(self,species_tree):
        
        pass

    def optimize_cost(self):
        pass

    
    def NNI(self):
        pass


