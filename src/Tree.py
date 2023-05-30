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
        




    def map(self):
        pass

    def optimize_cost(self):
        pass

    
    def NNI(self):
        pass


