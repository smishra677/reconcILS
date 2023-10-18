import utils.Tree as Tree
import argparse
import copy
import re
import os

class readWrite:
    def read_trees(self,i,folder):
        gene_tre= open(folder+'/rep_'+str(i)+'.tre')
        tr =gene_tre.read().strip().split('\n')
        gene_tre.close()
        return str(tr[0])

    def write_trees(self,i,folder,tree):
        with open(folder+'/rep_1_'+str(i)+'.tre', 'w') as f:
            f.write(tree)
        

    def read_log(self,flag,i,dic,folder):
        o= open('./'+folder+'/rep_'+str(i)+'.log').read()
        dic['Process']+=[flag]
        dic['Replicate']+=[i]
        


        for i in o.split('\n'):
            val= (i.split(':'))
            if val[0]=='Total duplications':
                dic['Duplication']+=[int(val[1])]
            elif val[0] in ['Rasmussen and Kellis hemiplasy']:
                dic['RHemiplasy']+=[int(val[1])]
            elif val[0] in ['All NNI']:
                dic['NNI']+=[int(val[1])]
            elif val[0] in ['All ILS (DLCPar)']:
                dic['DLCILS']+=[int(val[1])]
            elif val[0]=='Total losses':
                dic['Loss']+=[int(val[1])]
            elif val[0]=='Copy Number Hemiplasy':
                dic['Hemiplasy']+=[int(val[1])]
            else:
                continue
    
        return dic




    def Create_pd(self,flag,i,oo,dic):
        from collections import Counter
        o_list=[]

        for k in oo:
            if type(k)==list:
                o_list+=k
            else:
                o_list.append(k)

        o= dict(Counter(o_list))
                

        dic['DLCILS']+=[0]
        dic['RHemiplasy']+=[0]
        
        dic['Process']+=[flag]
        dic['Replicate']+=[i]
        dic['Duplication'].append(0)
        dic['NNI'].append(0)
        dic['Loss'].append(0)
        dic['Hemiplasy'].append(0)


        for i in o:
            if i in ['Duplication','NNI','Loss']:
                dic[i][-1]=o[i]



                
            
        

        return dic

    def Create_pd_ete3(self,flag,i,o,dic):
        
        dic['DLCILS']+=[0]
        dic['RHemiplasy']+=[0]
        
        dic['Process']+=[flag]
        dic['Replicate']+=[i]
        dic['Duplication'].append(0)
        dic['NNI'].append(0)
        dic['Loss'].append(0)
        dic['Hemiplasy'].append(0)


        for i in o:
            if i in ['D','L']:
                    if i =='D':
                        dic['Duplication'][-1]=o[i]
                    elif i=='L':
                        dic['Loss'][-1]=o[i]
                    else:
                        continue

        return dic
                        
                
                
            


    def sp_event_ete3(self,root):
        li=[]

        for node in root.traverse(strategy="postorder"):
            if len(node.children) != 0:
            ##print(dir(root))
                
                li.append(node.evoltype)
            ##print(root.isLeaf),
        
        return li


    def sp_event(self,root,li):
    
        if root:
            self.sp_event(root.leftChild,li)
            self.sp_event(root.rightChild,li)
            if (root.evolve=='NNI'):
                if root.cost in [0,1]:
                    li.append(root.evolve)
                else:
                    li+= ['NNI' for i in range(root.cost)]
            else:
                if type(root.evolve)==list:
                    for i in root.evolve:
                        if i=='NNI':
                            li+=['NNI' for k in range(root.cost)]
                        else:
                            li.append(i)
                else:
                    li.append(root.evolve)
            ##print(root.isLeaf),
        return li



    def to_newick(self,tree):
        newick = ""
        newick = self.traverse(tree, newick)
        newick = f"{newick};"
        return newick

    def traverse(self,tree, newick):
        if tree.leftChild and not tree.rightChild:
            newick = f"(,{self.traverse(tree.leftChild, newick)}){tree.taxa if tree.isLeaf else ''}"
        elif not tree.leftChild and tree.rightChild:
            newick = f"({self.traverse(tree.rightChild, newick)},){tree.taxa if tree.isLeaf else ''}"
        elif tree.leftChild and tree.rightChild:
            newick = f"({self.traverse(tree.rightChild, newick)},{self.traverse(tree.leftChild, newick)}){tree.taxa if tree.isLeaf else ''}"
        elif not tree.leftChild and not tree.rightChild :
            newick = f"{tree.taxa if tree.isLeaf else ''}"
        else:
            pass
        return newick




    def parse(self,newick):

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
                    new_tre= Tree.Tree()
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

        tre= Tree.Tree()
        val =recurse(tre)

        return val[-1]
