import os
import sys
sys.path.append("../../")
sys.path.append("../")
from reconcILS import *
from utils_reconcILS import *
import pandas as pd     
from collections import Counter




dic={'Process':[],'Replicate':[],'Gene_tree':[],'Species_Tree':[],'Duplication':[],'NNI':[],'Loss':[]}
    


#Species Tree String, alternatively you can read it as string from a file
sp_string='(A,(B,C));'






#Gene tree folder if you are using it to read the the gene tree 
gene_folder='./zero_loss'

#Change the cost of each events. 
Duplication_cost=1.1
ILS_cost=1.0
Loss_cost=1.0





#Gene tree as list of strings
list_gene_trees= ['((A,(B,C)),(B,(B,C)))']

#looping over all the gene trees
for i in range(len(list_gene_trees)):
        reconcILS =reconcils()
        red= readWrite.readWrite()
        reconcILS.D_cost=Duplication_cost
        reconcILS.I_cost=ILS_cost
        reconcILS.L_cost=Loss_cost
        

        
        #DATA_PATH = "./"
        #file_path = os.path.join(DATA_PATH)

   
        
        #Provide the gene tree as string. In this case we are getting gene trees from the list.
        gene_tree= str(list_gene_trees[i])
        


        sp_string = sp_string.replace('e-', '0')
        gene_tree = gene_tree.replace('e-', '0')



        tr= red.parse_bio(gene_tree)
        sp=red.parse_bio(sp_string)

        species_names,flag1,flag2=reconcILS.get_species(sp)


        if flag1:
            print('Error: please make sure all species are named with more than 2> or less than 2<= characters')
            exit()

        
     

        if not flag2:
            map_ = reconcILS.generate_letter_combinations()


            mapping = {name: next(map_) for name in species_names}


            resolve_sp_string=reconcILS.generate_species_codes(mapping,red.to_newick(sp))

            resolved_gene_tree=reconcILS.generate_species_codes(mapping,red.to_newick(tr))

            inverse_mapping = {value: key for key, value in mapping.items()}



            file_path_name_dict = './'+parser.output[:-4]+'species_acr.json'
            file_path_inverse_name_dict = './'+parser.output[:-4]+'acr_species.json'

            with open(file_path_name_dict, 'w') as file:
                json.dump(mapping, file)

            with open(file_path_inverse_name_dict, 'w') as file:
                json.dump(inverse_mapping, file)
        else:
            resolved_gene_tree=red.to_newick(tr)
            resolve_sp_string =red.to_newick(sp)


        sp=red.parse(resolve_sp_string)
        tr=red.parse(resolved_gene_tree)


        
        sp_copy= copy.deepcopy(sp)



        tr.order_gene(sp)
        tr.label_internal()
        sp.label_internal()
        tr.map_gene(sp)


        reconcILS.setCost(sp)
        sp.isRoot=True
        tr.isRoot=True
        sp_copy.isRoot=True
        

        reconcILS.gene_tree= copy.deepcopy(tr)
        re_w = readWrite.readWrite()
        species_edge_list=reconcILS.get_edges(sp)
        
        if reconcILS.introgression:
            red.write_introgression(sp)
        

        li= reconcILS.iterative_reconcILS(tr,sp,sp_copy,sp,[])
            
            
        
        #print('##############################################################################################')
        if len(li)==0:
            print('Error')
        else:   
            #print('##############################################################################################')
            li.reverse()
            dic_reconcILS={}
            for i in li:
                if type(i)==list:
                        if (repr(list(sorted(i[0]))),' to ',repr(list(sorted(i[1]))))  in dic_reconcILS.keys():
                            if i[3]==True:
                                dic_reconcILS[(repr(list(sorted(i[0]))),' to ',repr(list(sorted({i[1]}))))]+=[i[2]]
                            else:
                                dic_reconcILS[(repr(list(sorted(i[0]))),' to ',repr(list(sorted(i[1]))))]+=[i[2]]
                        else:   
                            if i[3]==True:
                                dic_reconcILS[(repr(list(sorted(i[0]))),' to ',repr(list(sorted({i[1]}))))]=[i[2]]
                            else:
                                dic_reconcILS[(repr(list(sorted(i[0]))),' to ',repr(list(sorted(i[1]))))]=[i[2]]

        #print('##############################################################################################')

        eve_rec= {'D':0,'I':0,'L':0}
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
                                    
      



                                    
                                    
        sp_small=red.parse(resolve_sp_string)
        sp_small.label_internal()


        reconcILS.edge_to_event(sp_small,dic_reconcILS,0) 

   



        if not flag2:
            dic['Gene_tree']+=[reconcILS.generate_code_gene(inverse_mapping,red.to_newick(reconcILS.gene_tree))]
            dic['Species_Tree']+=[reconcILS.generate_code_gene(inverse_mapping,red.to_newick_sp(sp_small))]

        else:
            dic['Gene_tree']+=[red.to_newick(reconcILS.gene_tree)]
            dic['Species_Tree']+=[red.to_newick_sp(sp_small)]


        dic['Replicate']+=[0]
        dic['Process']+=['reconcILS']  
        dic['Duplication']+=[eve_rec['D']]
        dic['NNI']+=[eve_rec['I']]
        dic['Loss']+=[eve_rec['L']]

        
            

        #dic= red.Create_pd('reconcILS',i,li,dic)


       


 

print(dic)
df = pd.DataFrame(dic)



#save the output as the file
df.to_csv(gene_folder+'_1_result.csv', index=False)



        


        


   
