
import Tree
import re

a= Tree.Tree()

tr= '((A,B),(D,C));'
sp ='((A,B),(D,C));'


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
 


def parse(newick):
    tokens = re.finditer(r"([^:;,()\s]*)(?:\s*:\s*([\d.]+)\s*)?([,);])|(\S)", newick+";")

    def recurse(tre,nextid = 0, parentid = -1): # one node
        thisid = nextid
        tre.id= nextid
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


printorder(tr)
print('###############################@22222222222222222222222222')
printorder(sp)