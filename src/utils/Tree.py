from itertools import count
import copy
import matplotlib.pyplot as plt
import pydot
import utils.idmaker_ as idmaker_




class Tree:
    def __init__(self):
        self.id= idmaker_.idmaker2().id
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
        self.event_list=[]
        self.visited=[]
        self.li=[]
        self.paralogy=[]
        self.NNI_=[]
        self.taxa_list=[]

    def clear_ref(self):
        if self:
            if self.leftChild:
                self.leftChild.clear_ref()
            if self.rightChild:
                self.rightChild.clear_ref()
            self.refTo=[]
            self.tag = None

    def reset(self):
        if self:
            self.event=None
            self.tag = None
            self.evolve=None
            self.taxa_list=[]
            self.NNI_=[]
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


    def visualize_binary_tree(self, graph):
        
        if self.leftChild:
            if len(str(self.leftChild.taxa))==0:
                label='None'
            else:
                label=str(self.leftChild.taxa)
            Node_left=pydot.Node(str(self.leftChild.id), label=label)
            node2 = graph.add_node(Node_left)
            if self.leftChild.evolve=='Loss':
                graph.add_edge(pydot.Edge(str(self.id),str(self.leftChild.id), color='blue'))
            else:
                graph.add_edge(pydot.Edge(str(self.id), str(self.leftChild.id), color='green'))
            self.leftChild.visualize_binary_tree(graph)
        if self.rightChild:
            if len(str(self.rightChild.taxa))==0:
                label='None'
            else:
                label=str(self.rightChild.taxa)
            Node_right=pydot.Node(str(self.rightChild.id), label=str(self.rightChild.taxa))
            node3 = graph.add_node(Node_right)
            if self.rightChild.evolve=='Loss':
                graph.add_edge(pydot.Edge(str(self.id),str(self.rightChild.id), color='blue'))
            else:
                graph.add_edge(pydot.Edge(str(self.id), str(self.rightChild.id), color='green'))
            self.rightChild.visualize_binary_tree(graph)

        else:
            return graph
        
        

        
        
        return graph


    def viz(self):
        graph = pydot.Dot()
        Node_root=pydot.Node(str(self.id), label=('None'))
        node1 = graph.add_node(Node_root)
        #graph.write_png("binary_tree.png")
        graph= self.visualize_binary_tree(graph)
        print(graph)
        graph.write_png("binary_tree.png")
       

    def find_all_edges_sp(root,edges,node_,taxa_,color_,LC_):

        node_.append(int(root.id))
        
        color_.append(root.evolve)
    
        LC_.append(root.taxa)
        
        
        if root.isLeaf:
            
            number_of_taxa=''.join(taxa_).count(root.taxa)
            taxa_.append(root.taxa+'_'+str(number_of_taxa))
            
        else:
            number_of_taxa=''.join(taxa_).count(root.evolve)
            taxa_.append(root.evolve+'_'+str(number_of_taxa+1))

        if root.rightChild:
            edges.append((root.id,root.rightChild.id))
            edges, node_,taxa_, color_,LC_=root.rightChild.find_all_edges(edges,node_,taxa_,color_,LC_)


        if root.leftChild:
            edges.append((root.id,root.leftChild.id))
            edges, node_,taxa_, color_,LC_ =root.leftChild.find_all_edges(edges,node_,taxa_,color_,LC_)

        return edges, node_,taxa_, color_,LC_




    def find_all_edges(root,edges,node_,taxa_,color_,LC_):
        node_.append(int(root.id))
        if type(root.evolve)==list:
            color_+=root.evolve
        else:
            if root.evolve==None:
                root.evolve =root.parent.evolve
                color_.append(root.parent.evolve)
            else:
                color_.append(root.evolve)   
        
        
        if root.isLeaf:
            LC_.append(set(root.taxa))
            number_of_taxa=''.join(taxa_).count(root.taxa)
            taxa_.append(root.taxa+'_'+str(number_of_taxa))
        else:
            LC_.append(root.taxa)
            if type(root.evolve)==list:
                for re in root.evolve:
                    number_of_taxa =''.join(taxa_).count(re)
                    taxa_.append(re+'_'+str(number_of_taxa+1))
            else:
                if root.evolve==None:
                    number_of_taxa=''.join(taxa_).count(root.parent.evolve)
                    taxa_.append(root.parent.evolve+'_'+str(number_of_taxa+1))
                else:
                    number_of_taxa=''.join(taxa_).count(root.evolve)
                    taxa_.append(root.evolve+'_'+str(number_of_taxa+1))
            #taxa_.append(root.evolve)
                 
      


        if root.rightChild:
            edges.append((root.id,root.rightChild.id))
            edges, node_,taxa_, color_,LC_=root.rightChild.find_all_edges(edges,node_,taxa_,color_,LC_)


        if root.leftChild:
            edges.append((root.id,root.leftChild.id))
            edges, node_,taxa_, color_,LC_ =root.leftChild.find_all_edges(edges,node_,taxa_,color_,LC_)


        #print(node_)
        #print(taxa_)
        #print(color_)

        return edges, node_,taxa_, color_,LC_
    
   
                
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
    

    def order_gene(self,species_tree):
    
        if self != None:
            if self.leftChild:
                self.leftChild.order_gene(species_tree)
            if  self.rightChild:
                self.rightChild.order_gene(species_tree)
            self.search_sp(species_tree),
    



    def label_internal(self):
        if self != None:
            if self.leftChild:
                self.leftChild.label_internal()
                if self.tag==None:
                    self.tag= [self.leftChild.tag]
                    self.taxa=set(self.leftChild.taxa)
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
                    self.taxa=self.taxa.union(self.rightChild.taxa)
                    if self.rightChild.isLeaf:
                        self.taxa_list.append(self.rightChild.taxa)
                        self.taxa_list.sort()

                    else:
                        self.rightChild.taxa_list.sort()
                        self.taxa_list.append(''.join(self.rightChild.taxa_list))  
                        self.taxa_list.sort()               

    


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
                    gene_tree.order_gene(new_recon_right)
                    gene_tree.label_internal()
                    new_recon_right.label_internal()
                    gene_tree.map_gene(new_recon_right)



                    gene_tree.reset()
                    gene_tree.order_gene(new_recon_left)
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
                    gene_tree.order_gene(self)
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


    

    def find_loss_sp(self,root):
        if self:
            #print('tag',self.to_newick())
            if self.leftChild and self.rightChild:
                if self.rightChild.isLeaf and len(self.rightChild.refTo)==0 and self.leftChild.isLeaf and len(self.leftChild.refTo)==0:
                    root.cost+=root.cost+1
                    return 1
            

            if self.leftChild:
                root.cost+=self.leftChild.find_loss_sp(root)
            if self.rightChild:
                root.cost+=self.rightChild.find_loss_sp(root)
            
            if len(self.refTo)==0:

                return 1
            else:
                return 0



                

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
                if self.evolve =='Loss':
                    self.cost=1
            else:
                self.cost = self.leftChild.cost+self.rightChild.cost
                if self.evolve in ['Loss']:
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


    def add_edges(self,graph,node):
        if node is not None:
            if node.leftChild is not None:
                graph.add_edge(node.taxa, node.leftChild.taxa)
                self.add_edges(graph, node.leftChild)
            if node.rightChild is not None:
                graph.add_edge(node.taxa, node.rightChild.taxa)
                self.add_edges(graph, node.rightChild)

    def plot_tree(root):
        graph = nx.DiGraph()
        root.add_edges(graph, root)
        pos = graphviz_layout(graph, prog='dot')
        nx.draw(graph, pos, with_labels=True, arrows=False)
        plt.show()




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
        tree2.order_gene(self)
        
        tree2.label_internal()
    
        self.label_internal()
        
        tree2.map_gene(self)
        node.label_internal()
        
        #val=self.find_total_cost(self)
        #print(self.find_cost(node,0))
        return self.find_cost(node,0)
        pass
    
    def write_events(self):
        from collections import Counter

        dic=dict(Counter(self.event_list))
        val= [str(k) + ':' + str(v) for k,v in dic.items()]

        return val if len(val)>0 else  ''

    def traverse(self, newick):
        if self.leftChild and not self.rightChild:
            newick = f"(,{self.leftChild.traverse(newick)}){self.taxa if self.isLeaf else self.write_events()}"
        elif not self.leftChild and self.rightChild:
            newick = f"({self.rightChild.traverse(newick)},){self.taxa if self.isLeaf else self.write_events()}"
        elif self.leftChild and self.rightChild:
            newick = f"({self.rightChild.traverse(newick)},{self.leftChild.traverse(newick)}){self.taxa if self.isLeaf else self.write_events()}"
        elif not self.leftChild and not self.rightChild :
            newick = f"{self.taxa if self.isLeaf else self.write_events()}"
        else:
            pass
        return newick




    def to_newick(self):
        newick = ""
        newick = self.traverse(newick)
        newick = f"{newick};"
        return newick



    def parse(self,newick):
        import re

        tokens = re.finditer(r"([^:;,()\s]*)(?:\s*:\s*([\d.]+)\s*)?([,);])|(\S)", newick+";")

        def recurse(tre,nextid = 0, parentid = -1): # one node
            thisid = nextid
            #tre.id= nextid
            children = []

            name, length, delim, ch = next(tokens).groups(0)
            tre.taxa= name
            if name!=0:
                tre.taxa= name.split('_')[0]
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

  
        

    def locate_copy(self,copy_tree,tree):
        if tree:

            if self.id==tree.id:

               

                if tree.parent.leftChild == tree:
                    
                    tree.parent.children.remove(tree.parent.leftChild)
                    tree.parent.leftChild =copy_tree
                    tree.leftChild =copy_tree.leftChild
                    tree.rightChild =copy_tree.rightChild
                    copy_tree.parent=tree.parent
                    tree.parent.children.append(copy_tree)

                    #print('left_child',tree.children)
                    return
                else:
                    
                    
                    tree.parent.children.remove(tree.parent.rightChild)
                    tree.parent.rightChild =copy_tree
                    tree.leftChild =copy_tree.leftChild
                    tree.rightChild =copy_tree.rightChild
                    copy_tree.parent=tree.parent
                    tree.parent.children.append(copy_tree)
                    #print('after1:',tree.parent.rightChild.to_newick())
                    #print('right_child',tree.children)
                    return
            else:
                self.locate_copy(copy_tree,tree.leftChild)
                self.locate_copy(copy_tree,tree.rightChild)                 
            


    def NNI1(self,gene_tree,node):

        geneTree_left =copy.deepcopy(gene_tree)
        geneTree_left.reset()

        geneTree_right =copy.deepcopy(gene_tree)
        geneTree_right.reset()

        copy_left = copy.deepcopy(self)
        copy_left.reset()

        add_left,add_right=None,None
        add_left = geneTree_left.find_address(node[0])

        add_right = geneTree_right.find_address(node[0])

        if node[2]=='Left':
            new_tree_left = Tree()
            new_tree_left.leftChild =node[1].leftChild
            new_tree_left.rightChild=node[0].rightChild
            

            new_tree_left.children =[new_tree_left.leftChild, new_tree_left.rightChild]
            new_tree_left.leftChild.parent=new_tree_left
            new_tree_left.rightChild.parent=new_tree_left
            

            add_left.rightChild= new_tree_left
            add_left.leftChild= node[1].rightChild


            new_tree_right = Tree()
            new_tree_right.leftChild =node[1].rightChild
            new_tree_right.rightChild=node[0].rightChild
            

            new_tree_right.children =[new_tree_right.leftChild, new_tree_right.rightChild]
            new_tree_right.leftChild.parent=new_tree_right
            new_tree_right.rightChild.parent=new_tree_right

            add_right.rightChild= new_tree_right
            add_right.leftChild= node[1].leftChild

        else:
            new_tree_left = Tree()
            new_tree_left.leftChild =node[1].leftChild
            new_tree_left.rightChild=node[0].leftChild
            

            new_tree_left.children =[new_tree_left.leftChild, new_tree_left.rightChild]
            new_tree_left.leftChild.parent=new_tree_left
            new_tree_left.rightChild.parent=new_tree_left
            

            add_left.leftChild= new_tree_left
            add_left.rightChild= node[1].rightChild


            new_tree_right = Tree()
            new_tree_right.leftChild =node[1].rightChild
            new_tree_right.rightChild=node[0].leftChild
            

            new_tree_right.children =[new_tree_right.leftChild, new_tree_right.rightChild]
            new_tree_right.leftChild.parent=new_tree_right
            new_tree_right.rightChild.parent=new_tree_right

            add_right.leftChild= new_tree_right
            add_right.rightChild= node[1].leftChild


        
        return [[self.parse(geneTree_left.to_newick()),self.parse(geneTree_left.to_newick()),'left'],[self.parse(geneTree_right.to_newick()),self.parse(geneTree_right.to_newick()),'right']]


    
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
            copy_left.children= [copy_left.rightChild, copy_left.leftChild]
            new_tree_left.parent=copy_left
            copy_right_child.parent=copy_left

            id_= copy_left_child.id
            new_tree_left.id= id_
            copy_left_child.id= new_tree_left.id

            new_tree_right = Tree()
            new_tree_right.leftChild=copy_right_child
            new_tree_right.rightChild=copy_right.rightChild
            new_tree_right.children = [new_tree_right.rightChild, new_tree_right.leftChild]
            new_tree_right.rightChild.parent=new_tree_right
            new_tree_right.leftChild.parent=new_tree_right




            copy_right.rightChild=new_tree_right
            copy_right.leftChild=copy_left_child
            copy_right.children= [copy_right.rightChild, copy_right.leftChild]
            copy_right.rightChild.parent=copy_right
            copy_right.leftChild.parent=copy_right


            id_= copy_left_child.id
            new_tree_right.id= id_
            copy_left_child.id= new_tree_right.id
        
        else:
            copy_right_child= copy.deepcopy(self.rightChild.rightChild)
            #print(copy_right_child.to_newick())
            copy_right_child.reset()
            copy_left_child= copy.deepcopy(self.rightChild.leftChild)
            #print(copy_left_child.to_newick())
            copy_left_child.reset()
            
            new_tree_left = Tree()
            new_tree_left.rightChild =copy_left_child
            new_tree_left.leftChild=copy_left.leftChild
            #print(new_tree_left.to_newick())
            new_tree_left.children =[new_tree_left.leftChild, new_tree_left.rightChild]
            new_tree_left.leftChild.parent=new_tree_left
            new_tree_left.rightChild.parent=new_tree_left


            copy_left.rightChild=new_tree_left
            copy_left.leftChild=copy_right_child
            copy_left.children= [copy_left.rightChild, copy_left.leftChild]
            copy_left.rightChild.parent=copy_left
            copy_left.leftChild.parent=copy_left

            id_= copy_left_child.id
            new_tree_left.id= id_
            copy_left_child.id= new_tree_left.id

            
            new_tree_right = Tree()
            new_tree_right.leftChild=copy_right_child
            new_tree_right.rightChild=copy_right.leftChild
            #print(new_tree_right.to_newick())
            new_tree_right.children= [new_tree_right.leftChild, new_tree_right.leftChild]
            new_tree_right.leftChild.parent=new_tree_right
            new_tree_right.leftChild.parent=new_tree_right

            copy_right.rightChild=new_tree_right
            copy_right.leftChild=copy_left_child
            copy_right.children= [copy_right.rightChild, copy_right.leftChild]
            copy_right.rightChild.parent=copy_right
            copy_right.leftChild.parent=copy_right        

            id_= copy_left_child.id
            new_tree_right.id= id_
            copy_left_child.id= new_tree_right.id

        if self.parent==None:
            #print('No_parent')
            #return [[self.parse(copy_left.to_newick()),self.parse(copy_left.to_newick()),'left'],[self.parse(copy_right.to_newick()),self.parse(copy_right.to_newick()),'right']]
            return [[copy_left,copy_left,'left'],[copy_right,copy_right,'right']]
        
        geneTree_left.label_internal()
        
        geneTree_right.label_internal()

        copy_left.label_internal()
        copy_right.label_internal()

        #geneTree_left.rightChild=copy_left

        #geneTree_right.leftChild=copy_right


        #print(geneTree_left)
        #print(geneTree_right)
        
        self.locate_copy(copy_left,geneTree_left)
        self.locate_copy(copy_right,geneTree_right)
        #print('Left:',geneTree_left.to_newick())
        #print('Right:',geneTree_right.to_newick())
        #return [[self.parse(copy_left.to_newick()),self.parse(geneTree_left.to_newick()),'left'],[self.parse(copy_right.to_newick()),self.parse(geneTree_right.to_newick()),'right']]
        return [[copy_left,geneTree_left,'left'],[copy_right,geneTree_right,'right']]