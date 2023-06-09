
import Tree
import re

a= Tree.Tree()

tr= '((A,C),B);'
sp ='((A,B),C);'
#tr='((B,C),(D,A));'

def tag(root):
 
    if root:
        tag(root.leftChild)
        print(root.taxa),
        print(root.event),
        if root.event !=None:
            print(root.event.refTo)
            print(root.event.taxa)
        #print(root.isLeaf),
        tag(root.rightChild)



def evolve(root):
 
    if root:
        evolve(root.leftChild)
        print(root.taxa),
        print(root.evolve),
        print(root.tag),
        #print(root.isLeaf),
        evolve(root.rightChild)


def sp_tag(root):
 
    if root:
        sp_tag(root.leftChild)
        print(root.taxa),
        print(root.refTo),
        #print(root.isLeaf),
        sp_tag(root.rightChild)

import copy

def sp_event(root):
 
    if root:
        sp_event(root.leftChild)
        print(root.taxa),
        print('event',root.evolve),
        print('cost',root.cost)
        #print(root.isLeaf),
        sp_event(root.rightChild)

def printInorder(root):
 
    if root:
        printInorder(root.leftChild)
        print(root.taxa),
        print(root.isLeaf),
        printInorder(root.rightChild)

def printorder(root):
 
    if root:
            print(root.tag),
            print(root.taxa),
            print(root.id),
            print(root.parent),
            print(root.refTo),
            print('###################'),
            printorder(root.leftChild)
            printorder(root.rightChild)
 

def to_newick(tree):
    newick = ""
    newick = traverse(tree, newick)
    newick = f"{newick};"
    return newick

def traverse(tree, newick):
    if tree.leftChild and not tree.rightChild:
        newick = f"(,{traverse(tree.leftChild, newick)}){tree.taxa if tree.isLeaf else ''}"
    elif not tree.leftChild and tree.rightChild:
        newick = f"({traverse(tree.rightChild, newick)},){tree.taxa if tree.isLeaf else ''}"
    elif tree.leftChild and tree.rightChild:
        newick = f"({traverse(tree.rightChild, newick)},{traverse(tree.leftChild, newick)}){tree.taxa if tree.isLeaf else ''}"
    elif not tree.leftChild and not tree.rightChild :
        newick = f"{tree.taxa if tree.isLeaf else ''}"
    else:
        pass
    return newick


def parse(newick):
    tokens = re.finditer(r"([^:;,()\s]*)(?:\s*:\s*([\d.]+)\s*)?([,);])|(\S)", newick+";")

    def recurse(tre,nextid = 0, parentid = -1): # one node
        thisid = nextid
        #tre.id= nextid
        children = []

        name, length, delim, ch = next(tokens).groups(0)
        tre.taxa= name
        if name!=0:
            tre.isLeaf= True
        if ch == "(":
            while ch in "(,":
                new_tre= Tree.Tree()
                node, ch, nextid,tre1 = recurse(new_tre,nextid+1, thisid)
                children.append(node)
                tre.children.append(tre1)

            tre.leftChild= tre.children[0]
            tre.children[0].parent =tre
            tre.children[1].parent =tre


            tre.rightChild=tre.children[1]
            name, length, delim, ch = next(tokens).groups(0)
        return {"id": thisid, "name": name, "length": float(length) if length else None, 
                "parentid": parentid, "children": children}, delim, nextid,tre

    tre= Tree.Tree()
    val =recurse(tre)

    return val[-1]


def find_parent_child(root,child):

    if len(root.refTo)>1:
            for tre in root.refTo:
                for tree1 in root.refTo:
                    if (tree1 in tre.children):
                        if tree1==tre.leftChild:
                            child.append([tre,tree1,'Left'])
                        else:
                            child.append([tre,tree1,'Right'])

    return child
        

            
def parent_child(root,child):
    if root:
        child= find_parent_child(root,child)
        child=find_parent_child(root.leftChild,child)
        child= find_parent_child(root.rightChild,child)
    return child

tr=parse(tr)
sp=parse(sp)
sp_copy= copy.deepcopy(sp)
sp_copy.reset()

tr.printorder_gene(sp)

tr.label_internal()
sp.label_internal()



tr.map_gene(sp)

new_gene_tree =copy.deepcopy(tr)
new_gene_tree.reset()
recon=Tree.Tree()
recon= copy.deepcopy(sp)
recon.tag_species(new_gene_tree,tr)
print(to_newick(recon))
if recon.split_list!=None:
    for i in recon.split_list:
        tr.map_recon(i)
else:
    tr.map_recon(recon)

print(to_newick(recon))

recon.clean_up()
recon.total_cost_()
sp_event(recon)

print(sp.find_cost(tr,0))

print('Gene_tree',to_newick(tr))





def ILS(gene_tree,tr,sp_copy,cost):
    list_tree=[]
    child=[]
    
   
    child= parent_child(tr,child)

    if len(child)==0 or cost==0:
        return gene_tree,cost
    else:
        child =[child[0]]
        for ch1 in child:
                ch = copy.deepcopy(ch1[0])
                ch.reset()

                geneTree =copy.deepcopy(gene_tree)
                geneTree.reset()
                list_tree= ch1[0].NNI(geneTree,ch1[2])
                best_cost=1000
                imporvement=False
                new_topo=copy.deepcopy(geneTree)
                for i in list_tree:
                    i[1].reset()
                    i[0].reset()
                    cop= copy.deepcopy(sp_copy)
                    cop.reset()
                        
                        
                    print('top_0',to_newick(i[0]))
                    print('topo_1',to_newick(i[1]))
                    print('##############################')
                    
                    new_cost =cop.optimize_cost(i[0],i[1])
   
                    #print(cost)
                    #print(new_cost,best_cost)
                    if best_cost>new_cost and cost>0:
                        best_cost=new_cost
                        new_topo=copy.deepcopy(i[1])
                        #print('new_topo',to_newick(new_topo))
                        imporvement=True
                
                if imporvement:
                        #print('new_topo1',to_newick(new_topo))

                        cost=cost-1
                        if cost==0:
                            return new_topo,cost
                        else:
                            new_sp = copy.deepcopy(sp_copy)
                            new_sp.reset()
                            new_topo.reset()
                            new_topo.printorder_gene(new_sp)
                        
                            new_topo.label_internal()
                    
                            new_sp.label_internal()
                        
                            new_topo.map_gene(new_sp)
                            return ILS(new_topo,new_sp,sp_copy,cost)

                        #print(cop.find_cost(tr))
                else:
                    return new_topo,cost
        

initial_cost=1
#sp.find_cost(tr,0)
new_topo,cost =(ILS(tr,sp,sp_copy,initial_cost))
print(to_newick(new_topo))
print(initial_cost-cost)




#tr.map_recon(recon)
print(to_newick(recon))
new_tr= copy.deepcopy(tr)
new_tr.reset()
print(new_tr.optimize_cost(recon,recon))
#recon.tag_loss()

print(to_newick(recon))