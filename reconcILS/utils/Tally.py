import utils.Tree as Tree
from collections import Counter
import pandas as pd
import igraph as ig
import matplotlib.pyplot as plt
    
class Tally:
    def __init__(self):
        self.edges=None
        self.node_=None
        self.taxa_=None
        self.color_=None
        self.LC_=None
        self.list_Duplication=None
        self.list_NNI=None
        

    def Duplication_NNI(self,sp,lis_paralogy,lis_NNI):
        if sp:
            list_paralogy,lis_NNI=self.Duplication_NNI(sp.leftChild,lis_paralogy,lis_NNI)
            list_paralogy,lis_NNI=self.Duplication_NNI(sp.rightChild,lis_paralogy,lis_NNI)
            if len(sp.paralogy)>0:
                lis_paralogy+=sp.paralogy
            if len(sp.NNI_)>0:
                lis_NNI+=sp.NNI_[0]
            return list_paralogy,lis_NNI

        return lis_paralogy,lis_NNI
    
    def get_Duplication_NNI(self,sp):
        self.edges, self.node_,self.taxa_, self.color_,self.LC_=sp.find_all_edges([],[],[],[],[])
        self.list_Duplication, self.list_NNI= self.Duplication_NNI(sp,[],[])


    def make_table(self,dic):
        print('##############################################################################')
            
        for i in self.list_Duplication:
            if(i[0]  in dic.keys() and i[1] in dic.keys()):
                print('#Duplication on edge between:',dic[i[0]],'---->',dic[i[1]],'#')
            else:
                for k in dic.keys():
                    print('#Duplication on edge between:',dic[i[0]],'---->',dic[i[1]],'#')
                    
                
        print('##############################################################################')
            
        for i in self.list_NNI:
            if(i[0]  in dic.keys() and i[1] in dic.keys()):
                print('#ILS on edge between:',dic[i[0]],'------->',dic[i[1]],'#')
            else:
                for k in dic.keys():
                    print('#ILS on edge between:',dic[i[0]],'------->',dic[i[1]],'#')
                        
        print('##############################################################################')
              
    def tally_NNI(self,tr,new_topo,orientation):
        if orientation=='Left':
            if new_topo.leftChild.taxa ==None:
                tr.NNI_+=[(tr.id,tr.id)]
                #tr.NNI_id_list+=[(tr,tr)]
            else:
                if len(new_topo.leftChild.taxa)==1:
                    if tr.leftChild.rightChild:            
                        tr.NNI_+=[(tr.id,tr.leftChild.id)]
                        #tr.NNI_id_list+=[(tr,tr.leftChild)]
                    else:
                        tr.NNI_+=[(tr.id,tr.rightChild.id)]
                        #tr.NNI_id_list+=[(tr,tr.rightChild)]
                else:
                     tr.NNI_+=[(tr.id,tr.id)]
                     #tr.NNI_id_list+=[(tr,tr)]
        else:
            if new_topo.rightChild.taxa==None:
                tr.NNI_+=[(tr.id,tr.id)]
                #tr.NNI_id_list+=[(tr,tr)]
            else:
                if  len(new_topo.rightChild.taxa)==1:
                    if tr.rightChild.leftChild:
                        tr.NNI_+=[(tr.id,tr.rightChild.id)]
                        #tr.NNI_id_list+=[(tr,tr.rightChild)]
                    else:
                        
                        tr.NNI_+=[(tr.id,tr.leftChild.id)]
                        #tr.NNI_id_list+=[(tr,tr.leftChild)]
                else:
                    tr.NNI_+=[(tr.id,tr.id)]
                    #tr.NNI_id_list+=[(tr,tr)]
                    



    def search_per(self,sp,recon):
        if sp:
            if sp.taxa==recon.taxa:
                sp.li.append(recon.evolve)
                if sp.isLeaf==None:
                    if sp.evolve==None:
                        sp.evolve=recon.evolve
                else:
                    sp.evolve='Speciation'
                
            self.search_per(sp.leftChild,recon)
            self.search_per(sp.rightChild,recon)

    def per(self,sp,recon):
        if recon:
            self.search_per(sp,recon)

            self.per(sp,recon.leftChild)
            self.per(sp,recon.rightChild)
        



    def make_graph(self,sp,sp_string,gene_tree):
        sp_1=sp.parse(sp_string)
        g_tree= sp_1.parse(gene_tree)
        g_tree.label_internal()
        sp.label_internal()

        #self.per(g_tree,sp)


        #sp.map_gene(sp_1)
        #sp= g_tree
        self.get_Duplication_NNI(sp)

        g = ig.Graph(edges=self.edges)
        new_dic_LC= dict(zip(self.node_,self.LC_))
        Keys = list(new_dic_LC.keys())
        Keys.sort()
        sorted_dict_LC = {i: new_dic_LC[i] for i in Keys}
        to_delete = [i for i in range(0,max(self.node_)) if i  not in self.node_]
        g.delete_vertices(to_delete)
        new_dic= dict(zip(self.node_,self.color_))
        Keys = list(new_dic.keys())
        Keys.sort()
        sorted_dict = {i: new_dic[i] for i in Keys}
        
        color_=list(sorted_dict.values())



        color_dict = {"Duplication": "blue", "Loss": "pink",'NNI':"red","Speciation":"yellow",None:"yellow"}
        g.vs["color"] = [color_dict[evolve] for evolve in color_]
   
    
   
        new_dic= dict(zip(self.node_,self.taxa_))

        #self.make_table(sorted_dict_LC)
        Keys = list(new_dic.keys())
        Keys.sort()
        sorted_dict = {i: new_dic[i] for i in Keys}
    
        g.vs["label"] = list(sorted_dict.values())
        
        s=str(g)

        fig, ax = plt.subplots()

        layout = g.layout_reingold_tilford(root=[int(s.split('\n')[3].split('--')[0])])
        ig.plot(g, layout=layout, target=ax)
        
        plt.show()
