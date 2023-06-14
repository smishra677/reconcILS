from itertools import count
import copy

class Tree:
    new_id = count()
    def __init__(self):
        self.id= next(Tree.new_id)
        self.taxa=None
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


    def reset_tag(self):
        if self:
            self.tag = None
    
            if self.leftChild:
                self.leftChild.reset()
            
            if self.rightChild:
                self.rightChild.reset()





    def reset_evolve(self):
        if self:
            if self.evolve=='Loss':
                self.evolve = None
    
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

            self.printorder_species(root.leftChild)
            self.printorder_species(root.rightChild)
            if self.isLeaf==True and (self.taxa== root.taxa) and self.tag == None:
                self.tag= root
                root.refTo.append(self)


    def search_sp(self,species_tree):
        self.printorder_species(species_tree)
        pass
    

    def printorder_gene(self,species_tree):
    
        if self != None:
            if self.leftChild:
                self.leftChild.printorder_gene(species_tree)
            if  self.rightChild:
                self.rightChild.printorder_gene(species_tree)
            self.search_sp(species_tree),
    



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
            if len(self.taxa)>1:
                if set(self.taxa).issubset((root.taxa)) and self.event == None:
                    self.event= root
                    root.refTo.append(self)



                    

    def search_spMap(self,species_tree):
        self.map_species(species_tree)
        pass

    def printorder(self):
    
        if self:
                print(self.taxa),
                print(self.evolve),
                print('###################'),
                if self.leftChild:
                    self.leftChild.printorder()
                if self.rightChild:
                    self.rightChild.printorder()


    def map_gene(self,species_tree):
        if self != None:
            if self.leftChild:
                self.leftChild.map_gene(species_tree)
            if  self.rightChild:
                self.rightChild.map_gene(species_tree)
            self.search_spMap(species_tree),
    

        pass
    

    def label_duplication_gene(self,tr):
        if tr:
            if self.taxa==tr.taxa:
                tr.evolve='Duplciation'
            else:
                self.label_duplication_gene(tr.leftChild)
                self.label_duplication_gene(tr.rightChild)

    
                
    def label_duplication(self,tr):
            for tre in self.refTo:
                    for tree1 in self.refTo:
                        if (tree1 in tre.children):
                            tre.label_duplication_gene(tr)
                        

    def tag_species(self,gene_tree,tr):
        if self:
            if  len(self.refTo)==0:
                    if self.leftChild:
                        self.leftChild.tag_species(gene_tree,tr)

                    if self.rightChild:
                        self.rightChild.tag_species(gene_tree,tr)
            if self.isLeaf:
                if len(self.refTo)>1 and self.parent.evolve not in ['Duplication','Speciation']:

                    self.refTo=[]
                    new_recon_right=copy.deepcopy(self)
                    
                    new_recon_right.reset()
                    new_recon_left=copy.deepcopy(self)
                    
                    new_recon_left.reset()
         
                    gene_tree.reset()
                    gene_tree.printorder_gene(new_recon_right)
                    gene_tree.label_internal()
                    new_recon_right.label_internal()
                    gene_tree.map_gene(new_recon_right)



                    gene_tree.reset()
                    gene_tree.printorder_gene(new_recon_left)
                    gene_tree.label_internal()
                    new_recon_left.label_internal()
                    gene_tree.map_gene(new_recon_left)
                    
                    new_recon_left.refTo=[]
                    new_recon_right.refTo=[]
                    self.evolve='Duplication'

                    
                    if self.leftChild:
                        #self.leftChild.evolve='Speciation'
                        self.leftChild=new_recon_left
                    if self.rightChild:
                        #self.rightChild.evolve='Speciation'
                        self.rightChild=new_recon_right
                    self.split_list= [new_recon_left,new_recon_right]
            
                else:
                    
                    if self.leftChild:
                        self.leftChild.tag_species(gene_tree,tr)


                    if self.rightChild:
                        self.rightChild.tag_species(gene_tree,tr)

                
            else:
                
                if len(self.refTo)>1:
                    self.refTo=[]
                    new_recon_right=copy.deepcopy(self)
                    
                    new_recon_right.reset()
                    new_recon_left=copy.deepcopy(self)
                    
                    new_recon_left.reset()

                    if self.leftChild:
                        self.leftChild=new_recon_left
                    if self.rightChild:
                        self.rightChild=new_recon_right
                    self.split_list= [new_recon_left,new_recon_right]


                    self.reset()
                    gene_tree.printorder_gene(self)
                    self.label_internal()
                    gene_tree.label_internal()
                    gene_tree.map_gene(self)
                    gene_tree.reset()
                    self.refTo=[]

                    self.leftChild.refTo=[]
                    self.rightChild.refTo=[]
                    self.evolve='Duplication'
                    self.label_duplication(tr)

                    

                    new_recon_left.tag_species(gene_tree,tr)
                    new_recon_right.tag_species(gene_tree,tr)

                else:
                    

                    if self.leftChild:
                        self.leftChild.tag_species(gene_tree,tr)

                    if self.rightChild:
                        self.rightChild.tag_species(gene_tree,tr)
        pass

    

    def tag_loss(self,root):
        if root:
            if root.isLeaf and not set(root.taxa).issubset((self.to_tag)) and root.evolve==None:
                root.evolve='Loss'
            self.tag_loss(root.leftChild)
            self.tag_loss(root.rightChild)
        pass
    

    def deactivate(self):
        if self:
            if self.evolve== None:
                self.evolve='Speciation'
            
            if self.leftChild:
                self.leftChild.deactivate()
            if self.rightChild:
                self.rightChild.deactivate()



    def clean_up(self):
        if self:
            if self.isLeaf and self.evolve==None:
                self.evolve='Loss'
            if self.leftChild:
                self.leftChild.clean_up()
            if self.rightChild:
                self.rightChild.clean_up()

    def search_sp_loss_(self,root):
        if root:
            self.search_sp_loss_(root.leftChild)
            self.search_sp_loss_(root.rightChild)

            if set(self.taxa).issubset((root.taxa)) and self.evolve == None:
                    self.evolve=root
                    if self.to_tag!=None:
                        self.tag_loss(root)
                    root.deactivate()

            else:
                if root.isLeaf==True:
                    if root.taxa not in list(self.taxa):
                        if  self.evolve==None and root.evolve==None:
                            root.evolve='Loss'

                else:
                    if self.taxa==root.taxa and root.evolve==None:
                        if self.to_tag!=None:
                            self.tag_loss(root)
                            root.deactivate()


    

                

    def search_sp_loss(self,recon):
        self.search_sp_loss_(recon)
        pass

    def map_recon(self,recon):
    
        if self != None:
            if self.leftChild:
                self.leftChild.map_recon(recon)
            if  self.rightChild:
                self.rightChild.map_recon(recon)

            if self.isLeaf == None:
                if self.leftChild.evolve!=None and self.leftChild.isLeaf==None:
                    self.to_tag=self.taxa.difference(self.leftChild.taxa)
                if self.rightChild.evolve!=None and self.rightChild.isLeaf==None:
                    self.to_tag=self.taxa.difference(self.rightChild.taxa)


                self.search_sp_loss(recon)






    def tag_gene(self):
        if self != None:
            if self.leftChild:
                self.leftChild.tag_gene()

            if self.rightChild:
                self.rightChild.tag_gene()
            left =set()
            right=set()
            if self.leftChild:
                left=set(self.leftChild.taxa)

            if self.rightChild:
                right = set(self.rightChild.taxa)
            if len(left & right  )>1:
                self.evolve='D'
            else:
                self.evolve='S'


           
            
    

        pass


    
    def total_cost_(self):
        if self:
            if self.leftChild:
                self.leftChild.total_cost_()
            if self.rightChild:
                self.rightChild.total_cost_()

            if self.isLeaf:
                if self.evolve !='Speciation':
                    self.cost=1
            else:
                self.cost = self.leftChild.cost+self.rightChild.cost
                if self.evolve in ['Duplication','NNI']:
                    self.cost = self.cost +1
                


    def sum_cost(self):
        if self:
            if self.leftChild:
                self.leftChild.total_cost_()
            if self.rightChild:
                self.rightChild.total_cost_()

            if not self.isLeaf:
                self.cost = self.leftChild.cost+self.rightChild.cost
                
                







    def sp_tag(self):
    
        if self:
            if self.leftChild:
                self.leftChild.sp_tag()
            print(self.taxa),
            print(self.refTo),
            #print(root.isLeaf),
            if self.rightChild:
                self.rightChild.sp_tag()







    def find_loss(self,recon_tree):
        recon_tree.map_gene(self)
        recon_tree.tag_loss()
        pass

    def find_cost(self,node,val):
        if self:
                        
            if self.leftChild:
                self.leftChild.find_cost(node,val)
            if self.rightChild:
                self.rightChild.find_cost(node,val)
            
            if self.taxa==node.taxa and len(self.refTo)>0:
                val= len(self.refTo)
                return val
            else:
                return val if val<15 else 15


        pass



    def find_total_cost(self,node):
        if node:
            val= len(node.refTo)
            left_val= self.find_total_cost(node.leftChild)
            right_val= self.find_total_cost(node.rightChild)
            if left_val>1:
                val= val+ left_val
            if right_val>1:
                val= val + right_val
            return val
        else:
            return 0

        pass

    def optimize_cost(self,node,tree2):
        tree2.printorder_gene(self)
        
        tree2.label_internal()
    
        self.label_internal()
        
        tree2.map_gene(self)
        node.label_internal()
        
        #val=self.find_total_cost(self)
        #print(self.find_cost(node,0))
        return self.find_cost(node,0)
        pass



    def locate_copy(self,copy_tree,tree):
        if tree:
            if tree.leftChild:
                self.locate_copy(copy_tree,tree.leftChild)
            
            if tree.rightChild:
                self.locate_copy(copy_tree,tree.rightChild)
            if tree.taxa==self.taxa:
                if tree.parent==None:
                    return
                if tree.parent.leftChild == tree:
                    tree.parent.children.remove(tree.parent.leftChild)
                    tree.parent.leftChild =copy_tree
                    tree.parent.children.append(copy_tree)

                    #print('left_child',tree.children)
                    return
                else:
                    tree.parent.children.remove(tree.parent.rightChild)
                    tree.parent.rightChild =copy_tree
                    tree.parent.children.append(copy_tree)
                    #print('right_child',tree.children)
                    return                  

            
    
    def NNI(self,gene_tree,flag):

        geneTree_left =copy.deepcopy(gene_tree)
        geneTree_left.reset()



        geneTree_right =copy.deepcopy(gene_tree)
        geneTree_right.reset()

        copy_left = copy.deepcopy(self)
        copy_left.reset()
        copy_right = copy.deepcopy(self)
        copy_right.reset()
        if flag=='Left':
            copy_right_child= copy.deepcopy(self.leftChild.rightChild)
            copy_right_child.reset()
            copy_left_child= copy.deepcopy(self.leftChild.leftChild)
            copy_left_child.reset()


            new_tree_left = Tree()
            new_tree_left.leftChild =copy_left_child
            new_tree_left.rightChild=copy_left.rightChild
            new_tree_left.children =[new_tree_left.leftChild, new_tree_left.rightChild]
            new_tree_left.leftChild.parent=new_tree_left
            new_tree_left.rightChild.parent=new_tree_left

            copy_left.rightChild=new_tree_left
            copy_left.leftChild=copy_right_child
            copy_left.children = [copy_left.rightChild, copy_left.leftChild]
            new_tree_left.parent=copy_left
            copy_right_child.parent=copy_left

            new_tree_right = Tree()
            new_tree_right.leftChild=copy_right_child
            new_tree_right.rightChild=copy_right.rightChild
            new_tree_right.children = [new_tree_right.rightChild, new_tree_right.leftChild]
            new_tree_right.rightChild.parent=new_tree_right
            new_tree_right.leftChild.parent=new_tree_right




            copy_right.rightChild=new_tree_right
            copy_right.leftChild=copy_left_child
            copy_right.children+= [copy_right.rightChild, copy_right.leftChild]
            copy_right.rightChild.parent=copy_right
            copy_right.leftChild.parent=copy_right
        
        else:
            copy_right_child= copy.deepcopy(self.rightChild.rightChild)
            copy_right_child.reset()
            copy_left_child= copy.deepcopy(self.rightChild.leftChild)
            copy_left_child.reset()
            
            new_tree_left = Tree()
            new_tree_left.rightChild =copy_left_child
            new_tree_left.leftChild=copy_left.leftChild
            new_tree_left.children =[new_tree_left.leftChild, new_tree_left.rightChild]
            new_tree_left.leftChild.parent=new_tree_left
            new_tree_left.rightChild.parent=new_tree_left


            copy_left.rightChild=new_tree_left
            copy_left.leftChild=copy_right_child
            copy_left.children+= [copy_left.rightChild, copy_left.leftChild]
            copy_left.rightChild.parent=copy_left
            copy_left.leftChild.parent=copy_left


            new_tree_right = Tree()
            new_tree_right.leftChild=copy_right_child
            new_tree_right.rightChild=copy_right.leftChild
            new_tree_right.children= [new_tree_right.leftChild, new_tree_right.leftChild]
            new_tree_right.leftChild.parent=new_tree_right
            new_tree_right.leftChild.parent=new_tree_right

            copy_right.rightChild=new_tree_right
            copy_right.leftChild=copy_left_child
            copy_right.children+= [copy_right.rightChild, copy_right.leftChild]
            copy_right.rightChild.parent=copy_right
            copy_right.leftChild.parent=copy_right        


        if self.parent==None:
            return [[copy_left,copy_left],[copy_right,copy_right]]
        
        geneTree_left.label_internal()
        
        geneTree_right.label_internal()

        copy_left.label_internal()
        copy_right.label_internal()


        self.locate_copy(copy_left,geneTree_left)
        self.locate_copy(copy_right,geneTree_right)
        return [[copy_left,geneTree_left],[copy_right,geneTree_right]]

