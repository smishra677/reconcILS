import dendropy
from dendropy.simulate import treesim
import numpy as np
import random
import sys
import copy
import io
import argparse
import os

# utility functions
def get_descendent_edges(edge):
    """
    Returns a list of all edges that "share" a node with ``self`` Only do head node, rather than head and tail as in get_adjacent_edges.
    """
    he = [i for i in edge.head_node.incident_edges() if i is not edge]
    return he

def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

def get_leaves(edge):
    leaves = []
    for leaf in edge.head_node.leaf_iter():
        leaves.append(str(leaf.taxon).strip("'"))
    return(leaves)

def get_all_leaves(tree):
    leaves = []
    for leaf in tree.leaf_iter():
        taxon = str(leaf.taxon).split()[0].strip("'")
        leaves.append(taxon)
    return(leaves)

def flatten(lst):
    result = []
    for element in lst:
        if isinstance(element, list):
            result.extend(flatten(element))
        else:
            result.append(element)
    return result

def parse_args():
    parser = argparse.ArgumentParser(description="Simulate gene family trees.")
    parser.add_argument('--stree', type=str, help="Path to file with newick formatted species tree to use for simulation. Branch lengths should be in coalescent units.")
    parser.add_argument('--mu_par', type=float, help="Loss rate.")
    parser.add_argument('--lambda_par', type=float, help="Duplication rate.")
    parser.add_argument('--reps', type=int, help="Number of gene families to simulate.")
    parser.add_argument('--output', type=str, help="Folder for storing results.")
    args= parser.parse_args()
    return(args)

def write_trees(tree, outputdir, rep):
    tree.write(path = '%s/rep_%s.tre' % (outputdir, rep), schema="newick", suppress_rooting=True)

def write_all_trees(trees, outputdir, rep):
    outfilename = '%s/rep_%s_alltrees_unmutated.tre' % (outputdir, rep)
    with open(outfilename, 'w') as f:
        for tree in trees:
            thestring = tree.as_string(schema="newick", suppress_rooting=True)
            f.write(thestring)

def write_mutated(trees, outputdir, rep):
    outfilename = '%s/rep_%s_alltrees_mutated .tre' % (outputdir, rep)
    with open(outfilename, 'w') as f:
        for tree in range(len(trees)):
            thestring = trees[tree].as_string(schema="newick", suppress_rooting=True)
            f.write(thestring)


def write_log(total_dups, total_losses, outputdir, rep, len_trees, table, annotated, hemiplasy, rk_hemiplasy, ils, ils_dlcpar, ils_joining, ils_joining_dlcpar):
    output_log_name = '%s/rep_%s.log' % (outputdir, rep)
    output_log = open(output_log_name, 'w')
    output_log.write("Replicate: %s\n" % rep)
    output_log.write("Total duplications: %s\n" % str(total_dups))
    output_log.write("Total losses: %s\n" % str(total_losses))
    output_log.write("Copy Number Hemiplasy: %s\n" % str(hemiplasy))
    output_log.write("Rasmussen and Kellis hemiplasy: %s\n" % str(rk_hemiplasy))
    #output_log.write("ILS: %s\n" % str(ils))
    #output_log.write("ILS (DLCPar): %s\n" % str(ils_dlcpar))
    #output_log.write("ILS joining: %s\n" % str(ils_joining))
    #output_log.write("ILS joining (DLCPar): %s\n" % str(ils_joining_dlcpar))
    output_log.write("All ILS: %s\n" % str(ils+ils_joining))
    output_log.write("All ILS (DLCPar): %s\n" % str(ils_dlcpar+ils_joining_dlcpar))
    if len_trees == 0:
        output_log.write("All copies were lost.\n")
    output_log.close()
    
    output_table_name = '%s/rep_%s.tsv' % (outputdir, rep)
    output_table = open(output_table_name, 'w')
    output_table.write(table)
    output_table.close()

    output_sptree_name = '%s/rep_sptree_%s.nex' % (outputdir, rep)
    output_sptree = open(output_sptree_name, 'w')
    output_sptree.write(annotated)
    output_sptree.close()

def check_lists_equality(list1, list2):
    return set(list1) == set(list2)

def get_gene_tree(annotated_sp_tree):
       # set up machinery for simulating gene trees
    gene_to_species_map = dendropy.TaxonNamespaceMapping.create_contained_taxon_mapping(
        containing_taxon_namespace=annotated_sp_tree.taxon_namespace,
        num_contained=1)
    
    # simulate a gene tree for each copy
    gene_tree = treesim.contained_coalescent_tree(containing_tree=annotated_sp_tree, gene_to_containing_taxon_map=gene_to_species_map)
    gene_tree.calc_node_root_distances()
    gene_tree.calc_node_ages()

    return(gene_tree)

def add_duplication(event_time, sp_leaves, new_subtree):

    # iterate over the edges of the gene tree. If an edge is a) alive at the right time and b) contains a relevant descedent, add it to the options
    potential_edges = []
    for gtedge in new_subtree.preorder_edge_iter():
        #print('\nProcessing a gene tree edge')
        if gtedge.tail_node != None:
            relevant_time_frame = [gtedge.head_node.age, gtedge.head_node.age + gtedge.length]
            #print(relevant_time_frame)
            if event_time > relevant_time_frame[0] and event_time < relevant_time_frame[1]:
                gtleaves = []
                for gtleaf in gtedge.head_node.leaf_iter():
                    gtleaves.append(str(gtleaf.taxon).split(' ')[0].strip("'"))
                if len(intersection(sp_leaves, gtleaves)) > 0:
                    potential_edges.append(gtedge)
                    #print(gtleaves)

    #print(len(potential_edges))

    # sample an edge
    the_edge = random.sample(potential_edges, 1)

    # extract the subtree
    leaves_to_sample = get_leaves(the_edge[0])
    subtree = new_subtree.extract_tree_with_taxa_labels(leaves_to_sample)
    subtree.calc_node_ages()

    for subtree_node in subtree.preorder_node_iter():
        subtending_length = event_time - subtree_node.age
        break

    # truncate the branch to the duplication time'
    for subtree_edge in subtree.preorder_edge_iter():
        subtree_edge.length = subtending_length
        break

    subtree.annotations['age'] = float(event_time)
    return(subtree)
    
def lose_a_copy(event_time, sp_leaves, all_trees, sp_tree):
    # find the potential edges that exist in all available trees
    potential_edges = []
    edge_trees = []
    edge_tree_indices = []
    tree_count = -1
    #print('The event time: %s' % event_time)
    #print(sp_leaves)

    for available_subtree in all_trees:
        #print(available_subtree)
                        
        # increment tree count
        tree_count += 1
                        
        for gtedge in available_subtree.preorder_edge_iter():
            #print('\nProcessing a gene tree edge')
            relevant_time_frame = [gtedge.head_node.age, gtedge.head_node.age + gtedge.length]
            #print(relevant_time_frame)
            #print(relevant_time_frame)
            if event_time > relevant_time_frame[0] and event_time < relevant_time_frame[1]:
                gtleaves = []
                for gtleaf in gtedge.head_node.leaf_iter():
                    gtleaves.append(str(gtleaf.taxon).split(' ')[0].strip("'"))
                if len(intersection(sp_leaves, gtleaves)) > 0:
                    potential_edges.append(gtedge)
                    edge_trees.append(available_subtree)
                    edge_tree_indices.append(tree_count)
                    #print(gtleaves)
                    
    #print(len(potential_edges))

    # sample an edge
    #print(len(potential_edges))
    the_edge_index = random.randint(0, len(potential_edges) - 1)
    the_edge = potential_edges[the_edge_index]
    the_tree = edge_trees[the_edge_index]
    the_tree_index = edge_tree_indices[the_edge_index]

    # check for hemiplasy on this edge.
    hemiplasy_count, mutation_dc_count = check_loss_hemiplasy(the_edge, sp_tree)
                
    # remove any daughters of the selected edge
    taxa_to_remove = []
    for item in the_edge.head_node.leaf_iter():
        taxa_to_remove.append(item.taxon)           

    try:
        updated_tree = the_tree.extract_tree_without_taxa(taxa_to_remove)
        try:
            tree_age = float(str(the_tree.annotations['age']).split("=")[1].strip("'"))
            updated_tree.annotations['age'] = float(tree_age)
            updated_tree.calc_node_ages()
        except:
            updated_tree.calc_node_ages()

        # replace the tree
        all_trees[the_tree_index] = updated_tree

    except:
        all_trees.pop(the_tree_index)
        #print('Throwing a tree out.')
    return(all_trees, hemiplasy_count, mutation_dc_count)

def get_adjacent_copy_numbers(adjacent_edges, all_trees):
    """Figure out how many copies are in each edge using the gene trees."""
    #print('Getting adjacent edges!')
    N_adjacent_edges = []
    for edge in adjacent_edges:
        edge_leaves = get_leaves(edge)
        #print('Edge!:')
        #print(edge_leaves)
        # how many of our trees have overlap in the descendents
        overlapping_trees = 0
        for tree in all_trees:
            #print('check tree:')
            #print(tree)
            gt_leaves = get_all_leaves(tree)
            if len(intersection(edge_leaves, gt_leaves)) > 0:
                overlapping_trees+=1
                #print('We have a copy!')
        N_adjacent_edges.append(overlapping_trees)
    #print(N_adjacent_edges)
    return(N_adjacent_edges)

def check_loss_hemiplasy(edge, sp_tree):
    
    dc_count = 0 # count CNH
    dc_2_count = 0 # count mutation deep coalescences


    leaves = get_leaves(edge)
    leaves = [x.split()[0] for x in leaves]

    time_of_duplication = edge.head_node.age + edge.length
        
    found_in_sp_tree = False
    exact_branch = False # use this to track whether the duplication is placed on the branch of the gene tree matching the branch of the species tree.

    for spedge in sp_tree.preorder_edge_iter():
        spleaves = get_leaves(spedge)
        current_check = check_lists_equality(spleaves, leaves)
        if current_check == True:
            found_in_sp_tree = True
            
        # check whether the matching branch exists at the time of duplication.
        try:
            sp_relevant_time_frame = [spedge.head_node.age, spedge.head_node.age + spedge.length]
        except: 
            sp_relevant_time_frame = [spedge.head_node.age, np.inf]

        if current_check == True and time_of_duplication > sp_relevant_time_frame[0] and time_of_duplication < sp_relevant_time_frame[1]:
            exact_branch = True
        
    if found_in_sp_tree == False:
        dc_count += 1
    elif exact_branch == False:
        dc_2_count += 1
    


    return(dc_count, dc_2_count)

def check_hemiplasy(mutated_subtree, sp_tree):
    
    dc_count = 0 # count CNH
    dc_2_count = 0 # count mutation deep coalescences

    for edge in mutated_subtree.preorder_edge_iter():

        leaves = get_leaves(edge)
        leaves = [x.split()[0] for x in leaves]

        time_of_duplication = edge.head_node.age + edge.length
        
        found_in_sp_tree = False
        exact_branch = False # use this to track whether the duplication is placed on the branch of the gene tree matching the branch of the species tree.

        for spedge in sp_tree.preorder_edge_iter():
            spleaves = get_leaves(spedge)
            current_check = check_lists_equality(spleaves, leaves)
            if current_check == True:
                found_in_sp_tree = True
            
            # check whether the matching branch exists at the time of duplication.
            try:
                sp_relevant_time_frame = [spedge.head_node.age, spedge.head_node.age + spedge.length]
            except: 
                sp_relevant_time_frame = [spedge.head_node.age, np.inf]

            if current_check == True and time_of_duplication > sp_relevant_time_frame[0] and time_of_duplication < sp_relevant_time_frame[1]:
                exact_branch = True
        
        if found_in_sp_tree == False:
            dc_count += 1
        elif exact_branch == False:
            dc_2_count += 1
    
        break



    return(dc_count, dc_2_count)

def tree_differences(sp_tree, gene_tree):
    ils_count = 0
    ils_dlcpar = 0
    for gtedge in gene_tree.postorder_edge_iter():
        match = False
        gt_relevant_leaves = get_leaves(gtedge)
        gt_relevant_leaves = [x.split()[0] for x in gt_relevant_leaves]
        for edge in sp_tree.postorder_edge_iter():
            relevant_leaves = get_leaves(edge)
            if check_lists_equality(relevant_leaves, gt_relevant_leaves) == True:
                match = True
        if match == False:
            ils_count += 1
    if ils_count > 0:
        ils_dlcpar = 1
    
    return(ils_count, ils_dlcpar)

def count_ils_daughter(gene_tree, sp_tree):
    # prune species tree
    leaves_to_keep = []
    for edge in gene_tree.postorder_edge_iter():
        leaves_to_keep.append(get_leaves(edge))
    leaves_to_keep = [item for row in leaves_to_keep for item in row]
    leaves_to_keep = set(leaves_to_keep)
    leaves_to_keep = list(leaves_to_keep)
    leaves_to_keep = [x.split()[0] for x in leaves_to_keep]

    species_subtree = sp_tree.extract_tree_with_taxa_labels(leaves_to_keep)
    species_subtree.calc_node_ages()

    ils_count = tree_differences(sp_tree=species_subtree, gene_tree=gene_tree)

    return(ils_count)           
                
 
def birth_death(sp_tree, lambda_par, mu_par):

    # intialize counters
    total_dups = 0
    total_losses = 0
    ils = 0 # our ils counter
    ils_dlcpar = 0 # dlcpar ils counter
    hemiplasy = 0 # CNH
    rk_hemiplasy = 0 # dc leading to mutation placement
    

    # sample the parent tree from the MSC
    parent_tree = get_gene_tree(sp_tree)
    #print("\n")
    #print("Species tree: ", sp_tree)
    #print("Gene tree: ", parent_tree)
    current_ils, current_ils_dlcpar = tree_differences(gene_tree=parent_tree, sp_tree=sp_tree)
    ils+=current_ils
    ils_dlcpar+=current_ils_dlcpar
    #print('This is my custom count', ils, ils_dlcpar)

    # set up list for storing trees, and add the parent tree
    all_trees = []
    unmodified_trees = []
    all_trees.append(parent_tree)
    unmodified_trees.append(parent_tree)

    # iterate over the edges of the species tree to add duplications and losses.
    for edge in sp_tree.preorder_edge_iter(): 

        # get leaves
        sp_leaves = get_leaves(edge)
        #print('We are on this edge: ')
        #print(sp_leaves)

        if edge.length == None: # skip root edge
            continue

        # add annotation of base copy number to the branch.
        if len(edge.annotations) == 0:
            edge.annotations['copies']='1'

        # get number of copies
        N=int(str(edge.annotations['copies']).split("=")[1].strip("'"))


        if N == 0:
            # stop because we start this branch with no copies.
            # update number of copies entering descendant branches. Use get_descendent_edges.
            adjacent_edges = get_descendent_edges(edge)

            # for each descendent edge, see how many subtrees are present and annotate accordingly
            adjacent_Ns = get_adjacent_copy_numbers(adjacent_edges, all_trees)

            # update numbers only for descendant nodes. 
            for edge in range(len(adjacent_edges)):
                adjacent_edges[edge].annotations['copies'] = adjacent_Ns[edge]

            continue


        #print("\nBeginning this branch with %s copies" % N)
        #print("This is a branch with length %s" % edge.length)
        
        # set tc to 0
        tc = 0.0

        while tc < edge.length and N > 0:

            # draw the time of the next event from an exponential distribution.
            current_scale_parameter =  1/((lambda_par + mu_par)*N)
            ti = np.random.exponential(scale = current_scale_parameter)

            # decide whether the event happens before the end of the branch
            if ti + tc > edge.length:
                N = N
                tc = edge.length


            else:
                tc = ti+tc

                # decide whether event is duplication or loss
                event_prob =  np.random.uniform(0,1)
                # calculate the event time (time from the base)
                event_time = (edge.length - tc) + edge.head_node.age

                if event_prob < (lambda_par / (lambda_par+mu_par)):
                    #print('duplication')

                    # draw a gene tree
                    new_subtree = get_gene_tree(sp_tree)
                    unmodified_trees.append(new_subtree)
                    #print(new_subtree)

                    # add the duplication
                    mutated_subtree = add_duplication(event_time, sp_leaves, new_subtree)
                    all_trees.append(mutated_subtree)
                    #print(mutated_subtree)

                    #check whether branches of the mutated subtree are present in the species tree.
                    hemiplasy_count,mutation_dc_count = check_hemiplasy(mutated_subtree, sp_tree)
                    hemiplasy += hemiplasy_count
                    rk_hemiplasy += mutation_dc_count
                    #print(mutated_subtree)
                    #print("This is my hemiplasy indicator: ", hemiplasy_count)
                    current_ils, current_ils_dlcpar = count_ils_daughter(gene_tree = mutated_subtree, sp_tree=sp_tree)
                    ils += current_ils
                    ils_dlcpar += current_ils_dlcpar
                    #print(sp_tree)
                    #print(mutated_subtree)
                    #print(current_ils, current_ils_dlcpar)
                    # increase number of copies
                    N = N + 1

                    # add info on duplication to edge annotations
                    edge.annotations['duplication_%s' % total_dups] = [(edge.length - tc) + edge.head_node.age]
                    total_dups += 1


                else:
                    #print('loss')

                    # decrease number of copies
                    N = N - 1

                    # decide what to lose and drop it
                    all_trees, hemiplasy_count, mutation_dc_count = lose_a_copy(event_time, sp_leaves, all_trees, sp_tree)
                    hemiplasy += hemiplasy_count
                    rk_hemiplasy += mutation_dc_count
                    #for tree in all_trees:
                    #    print(tree)


                    # add info on loss to edge annotations
                    edge.annotations['loss_%s' % total_losses] = [(edge.length - tc) + edge.head_node.age]

                    total_losses += 1

        # update number of copies entering descendant branches. Use get_descendent_edges.
        adjacent_edges = get_descendent_edges(edge)

        # for each descendent edge, see how many subtrees are present and annotate accordingly
        adjacent_Ns = get_adjacent_copy_numbers(adjacent_edges, all_trees)

        # update numbers only for descendant nodes. 
        for edge in range(len(adjacent_edges)):
            adjacent_edges[edge].annotations['copies'] = adjacent_Ns[edge]
    

    return(sp_tree, total_dups, total_losses, all_trees, unmodified_trees, hemiplasy, rk_hemiplasy, ils, ils_dlcpar)

def create_table(sp_tree, unmodified_trees, all_subtrees):

    # create annotated species tree
    count = 0


    for edge in sp_tree.preorder_edge_iter():
        edge.annotations.drop()
        edge.annotations['label'] = count
        count+=1
    newly_annotated = sp_tree.as_string(schema="nexus", suppress_annotations = False)


    



    # compare gene trees to the species tree
    all_branch_data = 'Branch\tParent_1'
    for item in range(len(unmodified_trees)-1):
        all_branch_data += '\t'
        all_branch_data += 'Copy_%s' % str(item+2)
    for edge in sp_tree.preorder_edge_iter():
        leaves = get_leaves(edge)
        if len(leaves) > 1 and len(leaves) < len(sp_tree.leaf_nodes()):
            #print("Look for edge:")
            #print(leaves)
            #print('check the gene trees')
            branch_data = '%s' % str(edge.annotations['label']).split("=")[1].strip("'")
            for genetree in unmodified_trees:
                #print(genetree)
                found_edge = False
                for gtedge in genetree.preorder_edge_iter():
                    if found_edge == False:
                        gtleaves = get_leaves(gtedge)
                        intersect_gtleaves = [x.split(' ')[0] for x in gtleaves]
                        if check_lists_equality(intersect_gtleaves, leaves):
                            found_edge = True
                branch_data += '\t'
                if found_edge:
                    branch_data += 'True'
                else:
                    branch_data += 'False'
            all_branch_data += '\n'
            all_branch_data += branch_data
    #print(len(unmodified_trees))

    return(all_branch_data, newly_annotated, unmodified_trees)
        
def check_ils_joining(the_coalesced_edge, subtreeleaves):

    ils = 0
    # joining leaves
    joining_leaves = [x.split()[0] for x in get_leaves(the_coalesced_edge)]
    subtreeleaves = [x.split()[0] for x in subtreeleaves]

    #print(subtreeleaves, joining_leaves)

    # check for match
    match = check_lists_equality(joining_leaves, subtreeleaves)
    if match == False:
        ils = 1

    return ils

def coalesce_subtrees(all_subtrees, annotated_sp_tree):

    ils_joining = 0
    ils_joining_dlcpar = 0

    # get the parent subtree
    parent_subtree = all_subtrees.pop(0)
    parent = parent_subtree
    #print('Here is the parent tree.')
    #print(parent)
    # list of available subtrees
    available_subtrees = []
    available_subtrees.append(parent_subtree)
    parent_height = parent.seed_node.distance_from_tip() + np.inf

    # sort subtrees by age
    ages = []
    for subtree in all_subtrees:
        age = float(str(subtree.annotations['age']).split("=")[1].strip("'"))
        ages.append(age)

    sorted_indices = sorted(range(len(ages)), key=lambda i: ages[i], reverse=True)

    sorted_subtrees = [all_subtrees[i] for i in sorted_indices]
    sorted_ages = [ages[i] for i in sorted_indices]



    # set the root edge in the species tree to infinity
    for edge in annotated_sp_tree.preorder_edge_iter():
        edge.length = np.inf
        break

    # set counter to keep track of copy number
    count = 1


    # iterate over the subtrees in decreasing order of age
    for thesubtree in range(len(sorted_subtrees)):
        #print('Coalesce this subtree')
        #print(sorted_subtrees[thesubtree])
        #print('To this parent!')
        #print(parent)

        # increment the copy number counter (first copy is labelled #2.)
        count += 1

        # grab the current subtree and the age of the subtree
        subtree = sorted_subtrees[thesubtree]
        age = sorted_ages[thesubtree]
        original_age = age

       # get the list of leaves in the subtree.
        subtree_leaves = subtree.leaf_nodes()
        subtree_leaves = [str(x.taxon).strip("'") for x in subtree_leaves]
     
        # Boolean to track whether we have managed to coalesce.
        coalesced = False

        # Use a while loop to find a place to coalesce the subtree.
        while coalesced == False:


            # Find the branch of the species tree on which coalescence should happen.
            species_tree_height = annotated_sp_tree.seed_node.distance_from_tip()

            # iterate over the edges, and check the timing of events against the age of the duplication and the leaves that should be present.
            for edge in annotated_sp_tree.levelorder_edge_iter():
                min_age = edge.head_node.distance_from_tip()
                try:
                    max_age = edge.tail_node.distance_from_tip()
                except:
                    max_age = np.inf
                if age >= min_age and age < max_age:
                    leaves = get_leaves(edge)
                    leaves_match = ['%s 1' % x for x in leaves]
                    if len(intersection(leaves_match, subtree_leaves))> 0:
                        # record the branch that is present at the correct time and has the correct leaves as the branch_to_coalesce.
                        branch_to_coalesce = [edge, leaves, min_age, max_age, max_age - age]
                        original_branch_leaves = leaves
            

            # Boolean to note whether we still need to look for the coalescence edge.
            find_edge_to_coal = True

           # while loop to find the edge of the parent gene tree that we should coalesce to.
            while find_edge_to_coal == True:

                # list to store edge heights (so we can know when first coalescent event happens)
                edge_heights = []
                possible_edges_to_coalesce = []

                # iterate over the edges, and check the timing of events against the age of the duplication and the leaves that should be present.
                for edge in parent.levelorder_edge_iter():
                    min_age = edge.head_node.distance_from_tip()
                    try:
                        max_age = edge.tail_node.distance_from_tip()
                    except:
                        max_age = parent_height

                    if age >= min_age and age < max_age:
                        leaves = get_leaves(edge)
                        intersect_leaves = [x.split(' ')[0] for x in leaves]
                        if len(intersection(intersect_leaves, branch_to_coalesce[1]))> 0:
                            possible_edges_to_coalesce.append([edge, leaves, min_age, max_age])
                            edge_heights.append(max_age)
                            
                # draw a coalescent time
                try:
                    combinations = len(possible_edges_to_coalesce)
                    tcoal = np.random.exponential(scale = 1/combinations)
                    my_edge_min= min(edge_heights)

                except:
                    find_edge_to_coal = False
                    tcoal = np.inf
                    my_edge_min = 0

                # do we coalesce before the shortest edge ends?
                if tcoal+age > my_edge_min and find_edge_to_coal == True:
                    branch_to_coalesce[-1] = branch_to_coalesce[3] - my_edge_min
                    if branch_to_coalesce[-1] <= 0:
                        # if we have failed and are now in a different edge of the species tree, we can break the parent gene tree loop and go get the new branch of the species tree.
                        find_edge_to_coal = False
                        tcoal = np.inf
                    else:
                        # otherwise, set the new age to the minimum edge height.
                        age = my_edge_min
                else:
                    # if we haven't failed, then we know which edge to coalesce to.
                    find_edge_to_coal = False


            # did coalescence happen in this branch?
            if tcoal < branch_to_coalesce[-1]:
                #print('we coalesced')

                # select an edge from available at random
                the_coalesced_edge_info = random.sample(possible_edges_to_coalesce, k=1)[0]
                the_coalesced_edge = the_coalesced_edge_info[0]

                # check for deep coalescence due to branch mismatch
                ils_joining_current = check_ils_joining(the_coalesced_edge, subtree_leaves)
                ils_joining += ils_joining_current
                ils_joining_dlcpar += ils_joining_current
    
                # update taxon namespace of new tree
                for leaf in subtree.leaf_node_iter():
                    number = sorted_indices[thesubtree] + 2
                    new_taxon = str(leaf.taxon.label).split()[0].strip("'") + ' '+ str(number)
                    leaf.taxon.label = new_taxon
                    leaf.taxon = dendropy.Taxon(new_taxon)
    
                # adjust branch length to account for any failures to coalesce to the initial parent branches
                for edge in subtree.preorder_edge_iter():
                    new_length = age - (edge.head_node.distance_from_tip()+edge.length) + edge.length
                    edge.length = new_length
                    break
                    

                # adjust branch length in subtree to include coalescence time
                for edge in subtree.preorder_edge_iter():
                    edge.length = edge.length + tcoal
                    total_edge = edge.length
                    total_height = edge.head_node.distance_from_tip() + edge.length
                    break

                # adjust branch length in parent tree
                for edge in parent.preorder_edge_iter():
                     if edge == the_coalesced_edge:
                        if edge.length == None:
                            edge.length = total_height - edge.head_node.distance_from_tip()
                            previous_edge = edge.length
                        elif edge.length < total_height - edge.head_node.distance_from_tip():
                            edge.length = total_height - edge.head_node.distance_from_tip()
                            previous_edge = edge.length   
                        else:
                            previous_edge = edge.length
                            edge.length = total_height - edge.head_node.distance_from_tip()

                        # branch lengths for daughter nodes
                        
                        this_new_length = edge.length
                        
                # get relevant nodes
                original_child = the_coalesced_edge.head_node
                original_parent_node = the_coalesced_edge.tail_node

                # create new subtree
                new_subtree = dendropy.Tree()
                for item in subtree.taxon_namespace:
                    new_subtree.taxon_namespace.add_taxon(item)
                for item in parent.taxon_namespace:
                    new_subtree.taxon_namespace.add_taxon(item)
                if previous_edge - this_new_length < 0:
                    sys.exit('noooo whyyy')
                cnode = new_subtree.seed_node.new_child(edge_length = previous_edge - this_new_length)
                cnode.add_child(subtree)
                cnode.add_child(original_child)

                # add new subtree to old subtree
                if the_coalesced_edge.head_node == parent.seed_node:
                    parent = new_subtree
                    #print('doing this')
                    #print(parent)
                else:
                    original_parent_node.remove_child(original_child)
                    original_parent_node.add_child(new_subtree)
                    #print('doing that')
                    #print(parent)

                coalesced = True
                
                # reformat tree
                output_stream = io.StringIO()
                sys.stdout = output_stream

                # convert the tree to a string and print it
                print(str(parent))
#
                # redirect stdout back to its original destination
                sys.stdout = sys.__stdout__
#
                # read the captured output from the IOStream object
                output = output_stream.getvalue()
                output = output + ';'

                # print the captured output
                parent = dendropy.Tree.get(data=output, schema="newick")
                parent.calc_node_ages()

    

            else:
                age = branch_to_coalesce[3]
                del branch_to_coalesce
    

    return(parent, ils_joining, ils_joining_dlcpar)
      
def main():
             
    # get arguments (species tree, lambda_par, mu_par, arbitrarily_long_root)
    args = parse_args()

    # create output folder
    if os.path.exists(args.output):
        sys.exit('The output folder already exists. Please remove, or change the output folder name and try again.')
    else:
        os.system('mkdir %s' % args.output)
        
    # specify species tree
    sp_tree_main = dendropy.Tree.get(path=args.stree, schema="newick")
    sp_tree_main.calc_node_root_distances()
    sp_tree_main.calc_node_ages()

    # set the arugments
    lambda_par = args.lambda_par
    mu_par = args.mu_par

    for i in range(args.reps):

        # copy of species tree
        sp_tree = sp_tree_main.clone()

        # perform top-down birth death
        annotated_sp_tree, dup_count, loss_count, all_trees, unmodified_trees, hemiplasy, rk_hemiplasy, ils, ils_dlcpar = birth_death(sp_tree, lambda_par, mu_par)

        # make table
        table, annotated_towrite, sorted_unmodified_subtrees = create_table(sp_tree.clone(), unmodified_trees, all_trees)
        write_all_trees(sorted_unmodified_subtrees, args.output, str(i))
        write_mutated(all_trees, args.output, str(i))


        if len(all_trees) == 0:
            print('replicate: %s' % i)
            ils_joining, ils_joining_dlcpar = [np.nan, np.nan]
            #print('All copies lost.')
            #write_log_alllost(dup_count, loss_count, args.output, str(i))


        else:
            print('replicate: %s' % i)

            # coalesce
            all_the_trees = [x.clone() for x in all_trees]
            simulated_gene_family_tree, ils_joining, ils_joining_dlcpar = coalesce_subtrees(all_the_trees, annotated_sp_tree)


            # print
            write_trees(simulated_gene_family_tree, args.output, str(i))
            del simulated_gene_family_tree

        write_log(dup_count, loss_count, args.output, str(i), len(all_trees), table, annotated_towrite, hemiplasy, rk_hemiplasy, ils, ils_dlcpar, ils_joining, ils_joining_dlcpar)

        # remove everything but the species tree
        del annotated_sp_tree
        del all_trees
        del sp_tree

if __name__ == "__main__":
    main()
