
import Tree
import re

a= Tree.Tree()

tr= '(((A,D),C),B);'
sp ='(((A,B),C),D);'
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
        print(root.event),
        print(root.refTo),
        #print(root.isLeaf),
        sp_event(root.rightChild)

def printInorder(root):
 
    if root:
        printInorder(root.leftChild)
        print(root.taxa),
        #print(root.isLeaf),
        printInorder(root.rightChild)

def printorder(root):
 
    if root:
            print(root.tag),
            print(root.taxa),
            print(root.id),
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
            tre.rightChild=tre.children[1]
            name, length, delim, ch = next(tokens).groups(0)
        return {"id": thisid, "name": name, "length": float(length) if length else None, 
                "parentid": parentid, "children": children}, delim, nextid,tre

    tre= Tree.Tree()
    val =recurse(tre)

    return val[-1]

tr=parse(tr)
sp=parse(sp)
sp_copy= copy.deepcopy(sp)
sp_copy.reset()

tr.printorder_gene(sp)

tr.label_internal()
sp.label_internal()



tr.map_gene(sp)


print(sp.find_cost(tr))

print('Gene_tree',to_newick(tr))



'''
def po(tr,sp_copy):
    list_tree=[]
    if len(tr.refTo)>1:
        child={}
        for tre in tr.refTo:
            child[tre]=[]
            for tree1 in tr.refTo:
                if (tree1 in tre.children):
                    child[tre].append(tree1)
        for ch1 in child:
            if len(child[ch1])>0:
                ch = copy.deepcopy(ch1)
                ch.reset()
                list_tree= ch.NNI('left')
                for i in list_tree:
                    i.reset()
                    cop= copy.deepcopy(sp_copy)
                    cop.reset()
                    
                    
                    print(to_newick(cop))
                    print(to_newick(i))
                    cop.optimize_cost(tr,i)
                    #print(cop.find_cost(tr))
        


po(sp,sp_copy)
'''