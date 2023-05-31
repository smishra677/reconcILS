
import Tree
import re

a= Tree.Tree()

tr= '((A,B),C);'
sp ='((A,C),B);'


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


def sp_tag(root):
 
    if root:
        sp_tag(root.leftChild)
        print(root.taxa),
        print(root.refTo),
        #print(root.isLeaf),
        sp_tag(root.rightChild)



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

'''
print('Inoder Traversal Gene tree')
print(printInorder(tr))
print('$$$$$$$$$$$$$$$')

print('Inoder Traversal Species tree')

print(printInorder(sp))
'''

print('###################################')
print('Mapping Leaf nodes from Gene to Species')
tr.printorder_gene(sp)



print('##################################')
print(tr.label_internal())
print(sp.label_internal())


#printorder(tr)
print('###############################@22222222222222222222222222')
#printorder(sp)


print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
tr.map_gene(sp)
#print(tag(tr))
#print(sp_tag(sp))

print('#####################55555555555555555555555555555555####')
recon=Tree.Tree()
recon=sp
recon.tag_species()

print(sp_event(sp))

print('#######################3')

printInorder(recon)
print(to_newick(recon))