from .utils_reconcILS import Tree
import argparse
import copy
import re
import os
import pandas as pd


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
        o= open(folder+'/rep_'+str(i)+'.log').read()
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

    def update_events(self,existing_events, new_events):
        for event in existing_events:
            for key in new_events:
                if key in event:
                    event[key] += 1

        return existing_events


    def write_introgression(self,tree):
        tree1 = copy.deepcopy(tree)
        branch= []
        events=[]
        stack = [tree1]
        while stack:
            curr=stack.pop()
            if curr.leftChild:
                if curr.leftChild.isLeaf==None:
                    branch.append([repr(list(sorted(curr.taxa))),' to ',repr(list(sorted(curr.leftChild.taxa)))])
                    lis_l =curr.NNI(curr,'Left')
                    lis_l[0][0].label_internal()
                    lis_l[1][0].label_internal()
                    introgression_l =sorted([sorted(list(lis_l[0][0].leftChild.taxa)),sorted(list(lis_l[0][0].rightChild.taxa))])
                    introgression_r =sorted([sorted(list(lis_l[1][0].leftChild.taxa)),sorted(list(lis_l[1][0].rightChild.taxa))])

                    events.append([{''.join(introgression_l[0])+'-'+''.join(introgression_l[1]):0,''.join(introgression_r[0])+'-'+''.join(introgression_r[1]):0}])
                    stack.append(curr.leftChild)
            if curr.rightChild:
                if curr.rightChild.isLeaf==None:
                    branch.append([repr(list(sorted(curr.taxa))),' to ',repr(list(sorted(curr.rightChild.taxa)))])
                    
                    lis_r =curr.NNI(curr,'Right')
                    lis_r[0][0].label_internal()
                    lis_r[1][0].label_internal()
                    introgression_l =sorted([sorted(list(lis_r[0][0].leftChild.taxa)),sorted(list(lis_r[0][0].rightChild.taxa))])
                    introgression_r =sorted([sorted(list(lis_r[1][0].leftChild.taxa)),sorted(list(lis_r[1][0].rightChild.taxa))])

                    events.append([{''.join(introgression_l[0])+'-'+''.join(introgression_l[1]):0,''.join(introgression_r[0])+'-'+''.join(introgression_r[1]):0}])
                    #events.append([lis_r[0][1],lis_r[1][1]])
                    stack.append(curr.rightChild)

        dic= {'branch':branch,'events':events}
        df = pd.DataFrame(dic)
        df.to_csv('./discordance.csv', index=False)

        

    def Create_pd(self,flag,i,oo,dic):
        from collections import Counter
        o_list=[]

        for k in oo:
            if type(k)==list:
                if k[2]=='I':

                    o_list+=['NNI']
                elif k[2]=='D':
                    o_list+=['Duplication']
                else:
                    o_list+=['Loss']
        

        o= dict(Counter(o_list))
                

        #dic['DLCILS']+=[0]
        #dic['RHemiplasy']+=[0]
        
        dic['Process']+=[flag]
        dic['Replicate']+=[i]
        dic['Duplication'].append(0)
        dic['NNI'].append(0)
        dic['Loss'].append(0)
        #dic['Hemiplasy'].append(0)


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
                        
                
                
    
    def search_ete3(self,root):

        stack=[]
        child_events=[]
        while stack:

            child_events += [node_child.evoltype for node_child in root.traverse(strategy="postorder") if len(node_child.children) >1]
                        
    def convert_1(self,df,df1):
        df['Loss'].append(sum(list(df1['loss'])))
        df['NNI'].append(sum(list(df1['coal'])))
        df['Duplication'].append(sum(list(df1['dup'])))

        
        return df

    def convert(self,df,df1):
        df['Loss'].append(sum(list(df1['loss'])))
        df['DLCILS'].append(sum(list(df1['coal'])))
        df['Duplication'].append(sum(list(df1['dup'])))
        df['Hemiplasy'].append(0)
        df['NNI'].append(0)
        df['RHemiplasy'].append(0)
        
        return df

    def put_it_in_sp(self,df,node_to_leaves,species_edge_list):

        events={idx:{'D':0,'I':0,'L':0} for idx in species_edge_list}
        for idx, row in df.iterrows():

            if int(row['parentid'])==-1:
                if row['nodeid'].isnumeric():
                    start= repr(list(sorted(node_to_leaves[str(row['nodeid'])])))
                    #events[(start,' to ',start)]= {'D':row['dup'],'I':row['coal'],'L':row['loss']}
                else:
                    start =repr(list({row['nodeid']}))
            
            else:
                start=repr(list(sorted(node_to_leaves[str(row['parentid'])])))
            end=''


            if row['nodeid'].isnumeric():
                end= repr(list(sorted(node_to_leaves[str(row['nodeid'])])))
            else:
                end =repr(list({row['nodeid']}))
            

            if (start,' to ',end) in events:
                events[(start,' to ',end)]= {'D':row['dup'],'I':row['coal'],'L':row['loss']}


        return events

    def label_reconcILS(self,li):
            li.reverse()
            dic_reconcILS={}
            for i in li:
                    #print(i)
                    #print(dic_reconcILS)
                    if type(i)==list:
                            if i[3]==True:
                                if (repr(list(sorted(i[0]))),' to ',repr(list(sorted({i[1]}))))  in dic_reconcILS.keys():
                                        #print(dic_reconcILS)
                                    dic_reconcILS[(repr(list(sorted(i[0]))),' to ',repr(list(sorted({i[1]}))))]+=[i[2]]
                                else:
                                    dic_reconcILS[(repr(list(sorted(i[0]))),' to ',repr(list(sorted({i[1]}))))]=[i[2]]
                                        
                            else:
                                if (repr(list(sorted(i[0]))),' to ',repr(list(sorted(i[1]))))  in dic_reconcILS.keys():
                                        #print(dic_reconcILS)
                                    dic_reconcILS[(repr(list(sorted(i[0]))),' to ',repr(list(sorted(i[1]))))]+=[i[2]]
                                else:
                                    dic_reconcILS[(repr(list(sorted(i[0]))),' to ',repr(list(sorted(i[1]))))]=[i[2]]
                
            return dic_reconcILS

    def label_ete3(self,recon_tree, fr):
        coun = 0
        coun1 = 0

        for node in recon_tree.traverse(strategy="preorder"):
            if len(node.children) >= 0:
                if hasattr(node, 'evoltype'):
                    if node.evoltype in ['D', 'L']:
                        if node.evoltype == 'D':
                            coun += 1

                        node_ = sorted(node.get_species())

                        if node.evoltype == 'L' and len(node.children) > 1:
                            if not hasattr(node.children[1], 'evoltype') or not hasattr(node.children[0], 'evoltype'):
                                continue
                            else:
                                if node.up:
                                    node_up = sorted(node.up.get_species())
                                    if node.evoltype == 'L' and node_up == node_:
                                        continue
                                    elif (repr(node_up), ' to ', repr(node_)) in fr.keys():
                                        fr[(repr(node_up), ' to ', repr(node_))] += [node.evoltype]
                                    else:
                                        fr[(repr(node_up), ' to ', repr(node_))] = [node.evoltype]
                                else:
                                    if (repr(node_), ' to ', repr(node_)) in fr.keys():
                                        fr[(repr(node_), ' to ', repr(node_))] += [node.evoltype]
                                    else:
                                        fr[(repr(node_), ' to ', repr(node_))] = [node.evoltype]

                                if hasattr(node.children[1], 'evoltype'):
                                    node.children[1].evoltype = 'S'
                                if hasattr(node.children[0], 'evoltype'):
                                    node.children[0].evoltype = 'S'

                        else:
                            if node.up:
                                node_up = sorted(node.up.get_species())
                                tem_node = node.up

                                while node_ == node_up:
                                    if tem_node.up:
                                        node_up = sorted(tem_node.up.get_species())
                                        tem_node = tem_node.up
                                    else:
                                        break

                                if node.evoltype == 'D':
                                    coun1 += 1

                                if (repr(node_up), ' to ', repr(node_)) in fr.keys():
                                    fr[(repr(node_up), ' to ', repr(node_))] += [node.evoltype]
                                else:
                                    fr[(repr(node_up), ' to ', repr(node_))] = [node.evoltype]
                            else:
                                if (repr(node_), ' to ', repr(node_)) in fr.keys():
                                    fr[(repr(node_), ' to ', repr(node_))] += [node.evoltype]
                                else:
                                    fr[(repr(node_), ' to ', repr(node_))] = [node.evoltype]
        return fr





    def sp_event_ete3(self,root,events):
        li=[]
        

        
        for ev in events:
            if ev.etype == "D": 
                #print(ev.inparalogs, ev.outparalogs)
                li.append(ev.etype)

        #print(li)
        for node in root.traverse(strategy="preorder"):
            
            if hasattr(node,'evoltype') :
                #print(node)
                #print(node.evoltype)
                #print('############################')
                if (node.evoltype)!='D':
                    node_up= node.up
                    #print(node)

                    #print(node_up)
                   
                    #print('$$$$$$$$$$$$$$$$$$$$$$$')
                    #print(child_events)
                    if node.evoltype=='L':


                        child_events = []
                        for node_child in node.traverse(strategy="preorder"):
                            if hasattr(node_child,'evoltype'):
                                if  node!=node_child:
                                    child_events.append(node_child.evoltype)
                            else:
                                child_events.append('Z')


                        #print(node,child_events)
                    
                        
                        if 'S' in child_events or 'D' in child_events or 'X' in child_events or 'Z' in child_events :
                            if node.evoltype=='D':
                                print('Wrong Marking X')
                                exit()
                            #print('set to X-------------------------------------------------------')
                            setattr(node, 'evoltype', 'X')
                        elif len(child_events)==0:
                            #print(node)
                            li.append(node.evoltype)

                        else: 
                            #print(True)

                            #print(node)
                            li.append(node.evoltype)
                            for node_child in node.traverse(strategy="preorder"):
                                if hasattr(node_child, 'evoltype'):
                                    if node_child.evoltype=='L'and node_child!=node:
                                        if node_child.evoltype=='D':
                                            print('wrong marking')
                                            exit()
                                        setattr(node_child, 'evoltype', 'X')
                                        #print('set to Y-------------------------------------------------------')

                            




            ##print(dir(root))
                
                    #li.append(node.evoltype)
            ##print(root.isLeaf),
        
        return li

    def find_leaves(self,node, df):
        if not node.isnumeric():
            return node


        children = df[df['parentid'] == int(node)]['nodeid'].tolist()

        if not children:
            return {node}
        else:
            result = set()
            for child in children:
                result.update(self.find_leaves(child, df))
            return result


    def Create_pd_ete3_1(self,flag,i,o,dic):
            dic['Process']+=[flag]
            dic['Replicate']+=[i]
            dic['Duplication'].append(0)
            dic['Loss'].append(0)
            dic['NNI'].append(0)

            

            for i in o:
                #print(i)
                if i in ['D','L']:
                        if i =='D':
                            dic['Duplication'][-1]=o[i]
                        elif i=='L':
                            dic['Loss'][-1]=o[i]
                        else:
                            continue

            return dic

    def find_replacement_for_edge(self,species_edge_list,fr):
            list_to_replace ={}
            for edge in fr:
                    if edge not in species_edge_list:
                            for edge_sp in species_edge_list:
                                    if edge[2] in edge_sp[2]:
                                            list_to_replace[edge]=edge_sp


            for edge in list_to_replace:
                    fr[list_to_replace[edge]]= fr[edge]
                    del fr[edge]
            
            return fr

    def put_it_in_sp_tree(self,species_edge_list,dic_reconcILS,fr):
            eve_rec= {'D':0,'I':0,'L':0}
            eve_ete= {'D':0,'I':0,'L':0}
            from collections import Counter
            for i in species_edge_list:
                    rec_in_dic= {'D':0,'I':0,'L':0}
                    if i in dic_reconcILS:
                            dic_in= dict(Counter(dic_reconcILS[i]))
                            for j in rec_in_dic:
                                    if j in dic_in:
                                            rec_in_dic[j]=dic_in[j]
                                            eve_rec[j]+=dic_in[j]
                            dic_reconcILS[i]=rec_in_dic
                    else:
                            dic_reconcILS[i]={'D':0,'I':0,'L':0}
                                                            
                                                    
                    ete_in_dic= {'D':0,'I':0,'L':0}

                    if i in fr:
                            dic_in= dict(Counter(fr[i]))
                            for j in ete_in_dic:
                                    if j in dic_in:
                                            ete_in_dic[j]=dic_in[j]
                                            eve_ete[j]+=dic_in[j]
                            fr[i]=ete_in_dic
                    else:
                            fr[i]={'D':0,'I':0,'L':0}
            return dic_reconcILS,fr,eve_rec,eve_ete
                    



    def sp_event_gene(self,root,li):
        from collections import Counter
    
        if root:
            self.sp_event_gene(root.leftChild,li)
            self.sp_event_gene(root.rightChild,li)
            li.append([root.taxa,dict(Counter(root.event_list[1]))])
            
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
    
    def label_lost_child(self,tree):
        if tree:
            tree.event_list.append([-2,['L']])
            #print(tree.taxa)
            #print(tree.isLeaf)
            #print(tree.leftChild)
            #print(tree.rightChild)
            #print(tree)
            print('---------------------------')

        
            self.label_lost_child(tree.leftChild)

    


            self.label_lost_child(tree.rightChild)





    def write_events(self,tree):
        from collections import Counter
        if tree.isLeaf:
            if len(tree.event_list)>0:
                ev=''
                tree.event_list.reverse()
                tree.event_list=self.summarize(tree.event_list)
                for j in tree.event_list:
                    if j[0] not in [-1,-2]:
                        dic=dict(Counter(j[1]))
                        val= [str(v)+str(k) for k,v in dic.items()]
                        if len(ev)>0:
                            ev+='+'+''.join(val)+'  '+j[0]
                        else:
                            ev+=''.join(val)+'  '+j[0]
                    else:
                        dic=dict(Counter(j[1]))
                        val= [str(k) for k,v in dic.items()]
                        if len(ev)>0:
                            ev+='+'+'+'.join(val)
                        else:
                            ev+='+'.join(val)

                
                return (ev+'  '+tree.numbered_taxa)
            else:
                return tree.numbered_taxa+ ''
            
            
        else:
            if len(tree.event_list)>0:
                ev=''
                tree.event_list.reverse()
                tree.event_list=self.summarize(tree.event_list)
                #print(tree.event_list)
                for j in tree.event_list:
                    if j[0] not in [-1,-2]:
                        dic=dict(Counter(j[1]))
                        val= [str(v)+str(k) for k,v in dic.items()]
                        
                        if len(ev)>0:
                            ev+='+'+''.join(val)+'  '+j[0]
                        else:
                            ev+=''.join(val)+'  '+j[0]
                    else:
                        dic=dict(Counter(j[1]))
                        val= [str(k) for k,v in dic.items()]
                        if len(ev)>0:
                            ev+='+'+'+'.join(val)
                        else:
                            ev+='+'.join(val)
  

  

                    

                return ev
            else:
                return ''
    


    def write_events_sp(self,tree):
        if tree.isLeaf:
            ev=''
            for i in tree.event_list:
                ev+=str(i[0]['D'])+'-'+str(i[0]['I'])+'-'+str(i[0]['L'])
                ev+='   '
            return ev+tree.taxa   
        else:
            ev=''
            for i in tree.event_list:
                ev+=str(i[0]['D'])+'-'+str(i[0]['I'])+'-'+str(i[0]['L'])
                ev+='   '
                
            return ev

    def summarize(self,li):
        changed_list = [li[0]]

        for i in li[1:]:
            if i[0] == changed_list[-1][0]:
                changed_list[-1][1].extend(i[1])
            else:
                changed_list.append(i)

        return changed_list

    # Got it From StackOverflow:
    #https://stackoverflow.com/questions/61117131/how-to-convert-a-binary-tree-to-a-newick-tree-using-python
    # https://stackoverflow.com/questions/61117131/how-to-convert-a-binary-tree-to-a-newick-tree-using-python

    def to_newick(self,tree):
        newick = ""
        newick = self.traverse(tree, newick)
        newick = f"{newick};"
        return newick
    

    # Got it From StackOverflow:
    #https://stackoverflow.com/questions/61117131/how-to-convert-a-binary-tree-to-a-newick-tree-using-python
    # https://stackoverflow.com/questions/61117131/how-to-convert-a-binary-tree-to-a-newick-tree-using-python

    def traverse(self,tree, newick):
        if tree.leftChild and not tree.rightChild:
            newick = f"(,{self.traverse(tree.leftChild, newick)}){self.write_events(tree)}"
        elif not tree.leftChild and tree.rightChild:
            newick = f"({self.traverse(tree.rightChild, newick)},){self.write_events(tree)}"
        elif tree.leftChild and tree.rightChild:
            newick = f"({self.traverse(tree.rightChild, newick)},{self.traverse(tree.leftChild, newick)}){self.write_events(tree)}"
        elif not tree.leftChild and not tree.rightChild :
            newick = f"{self.write_events(tree)}"
        else:
            pass
        return newick



    # Got it From StackOverflow:
    #https://stackoverflow.com/questions/61117131/how-to-convert-a-binary-tree-to-a-newick-tree-using-python
    # https://stackoverflow.com/questions/61117131/how-to-convert-a-binary-tree-to-a-newick-tree-using-python
    def traverse_sp(self, newick,tree):
        if tree.leftChild and not tree.rightChild:
            newick = f"(,{self.traverse_sp(newick,tree.leftChild)}){self.write_events_sp(tree)}"
        elif not tree.leftChild and tree.rightChild:
            newick = f"({self.traverse_sp(newick,tree.rightChild)},){self.write_events_sp(tree)}"
        elif tree.leftChild and tree.rightChild:
            newick = f"({self.traverse_sp(newick,tree.rightChild)},{self.traverse_sp(newick,tree.leftChild)}){self.write_events_sp(tree)}"
        elif not tree.leftChild and not tree.rightChild :
            newick = f"{self.write_events_sp(tree)}"
        else:
            pass
        return newick





    # Got it From StackOverflow:
    #https://stackoverflow.com/questions/61117131/how-to-convert-a-binary-tree-to-a-newick-tree-using-python
    # https://stackoverflow.com/questions/61117131/how-to-convert-a-binary-tree-to-a-newick-tree-using-python
    def to_newick_sp(self,tree):
        newick = ""
        newick = self.traverse_sp(newick,tree)
        newick = f"{newick};"
        return newick

    #Got the framework of the code from StackOver FLow 
    # https://stackoverflow.com/questions/51373300/how-to-convert-newick-tree-format-to-a-tree-like-hierarchical-object
    # https://stackoverflow.com/questions/51373300/how-to-convert-newick-tree-format-to-a-tree-like-hierarchical-object



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
                tre.numbered_taxa=name
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


   

    def parse_bio(self, newick):
        import random
        tokens = re.finditer(r"([^:;,()\s]*)(?:\s*:\s*([\d.]+)\s*)?([,);])|(\S)", newick + ";")
        random.seed(42)
        def recurse(tre, nextid=0, parentid=-1): 
            thisid = nextid
            children = []

            name, length, delim, ch = next(tokens).groups(0)
            tre.taxa = name
            if name:
                tre.taxa = ''.join(name.split('_')[:2])
                tre.numbered_taxa = name
                tre.isLeaf = True

            if ch == "(":
                while ch in "(,":
                    new_tre = Tree.Tree()
                    node, ch, nextid, tre1 = recurse(new_tre, nextid + 1, thisid)
                    children.append(node)
                    tre.children.append(tre1)

                if len(tre.children) > 2:
                    random.shuffle(tre.children)
                    
                    
                    tre.leftChild = tre.children[0]
                    tre.rightChild = tre.children[1]
                    tre.children[0].parent = tre
                    tre.children[1].parent = tre
                    current = tre.rightChild
                    for child in tre.children[2:]:
                        ne_t = current.deepcopy()
                        child.parent = current
                        current.rightChild = child
                        current.leftChild =ne_t
                        ne_t.parent=current
                        current.taxa=''
                        current.isLeaf=False
                        current = child

                elif len(tre.children) == 1:
                    tre.children =[]
                    tre =tre1
                elif len(tre.children) == 2:
                    tre.leftChild = tre.children[0]
                    tre.rightChild = tre.children[1]
                    tre.children[0].parent = tre
                    tre.children[1].parent = tre

                name, length, delim, ch = next(tokens).groups(0)

            return {"id": thisid, "name": name, "length": length if length else None,
                    "parentid": parentid, "children": children}, delim, nextid, tre

        tre = Tree.Tree()
        val = recurse(tre)

        return val[-1]
