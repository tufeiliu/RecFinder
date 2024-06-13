from ..utils import files, util, process_args, parallel_task_manager
import dendropy
from typing import List, Set, Dict, Tuple
import pandas as pd
import numpy as np
from tqdm import tqdm
import multiprocessing as mp
import warnings
import scipy as spy
import os
import threading
import traceback
try:
    import queue
except ImportError:
    import Queue as queue     


def get_child_parent_relationships(treefile: str, outgroup_mrca: List[str], mrca_of_interest: List[str]):

    '''INPUT: treefile - phylogenetic tree in newick format
              outgroup_mrca - list of species in outgroup (required to root the tree)
              mrca_of_interest - 
       RETURNS: list of child -> parent relationships within the tree
    '''
    # read in phylogenetic tree
    t = dendropy.Tree.get(file=open(treefile, 'r'), schema="newick")
    # identify the most recent common ancestor (mrca) of the list of species in the outgroup_mrca list
    outgroup_mrca = t.mrca(taxon_labels=outgroup_mrca)
    # get the name of the mrca node
    outgroup_label = outgroup_mrca.label
    # re-root the tree at the node edge of the outgroup
    t.reroot_at_edge(outgroup_mrca.edge, update_bipartitions=False)
    outgroup_mrca.label = outgroup_label
    
    #initiate the list of child, parents of each branch in the tree 
    child_parent_set = set()
    root = ''
    internal_node_count = 0
    taxon_count = 0
    
    # get the most recent common ancestor of the tree we want to count (note: there needs to be a different way of doing this)
    tree_mrca = t.mrca(taxon_labels=mrca_of_interest)
    root_node = tree_mrca.label
    #iterate through the nodes in the clade defined by mrca_of_interest
    for nd in tree_mrca.postorder_iter():

        
        if nd.label == tree_mrca.label:
            continue
        
        internal_node_count +=1
        
        #identify the root (has no child -> parent relationship)
        if nd.parent_node == None:
            root = nd.label
            continue

        # get child-parent relationships for internal nodes and leaves
        if nd.is_internal():  # Internal node
            child_label = nd.label
        else:  # Leaf node
            taxon_count += 1
            child_label = nd.taxon.label
        
        # Add the child-parent relationship
        if child_label is not None and nd.parent_node.label is not None:
            child_parent = (child_label, nd.parent_node.label)
            child_parent_set.add(child_parent)
        
        # #get child-parent relationships for internal nodes and leaves
        # elif nd.parent_node != None: #confirming a true child -> parent relationship
            
        #     #child is an internal node (use nd.label for child)
        #     if nd.label != None: #internal nodes have a node label (leaves seemingly do not)
        #         # child = nd.label
        #         child_parent = (nd.label, nd.parent_node.label)
        #     #child is a leaf (use nd.taxon.label for child)
        #     else:
        #         taxon_count += 1
        #         # child = nd.taxon.label #get child -> parent relationships for leaves
        #         child_parent = (nd.taxon.label, nd.parent_node.label)
        # # add the child, parent relationship for the branch 
        # child_parent_set.add(child_parent)
    
    # print(f"The number of nodes in the tree analysed is {internal_node_count}") # = (n-1), where n = number of species 
    # print(f"The number of leaves in the tree analysed is {taxon_count}") # number of species analysed
    return child_parent_set, root_node


def CountMutations(child_parent_set: Set[Tuple[str, str]], 
                   sequence_dict: Dict[str, str], 
                   residue_dict: Dict[str, str], 
                   residue_loc: int):
    '''
    INPUT: 
    RETURNS: site matrix (20x20) with all substitutions at a specific site
    '''

    residue_mut = np.zeros((20, 20), dtype=np.int64)
    #iterate over every branch in the tree
    for edge in child_parent_set:
        
        child = edge[0]
        parent = edge[1]
        child_ident = sequence_dict[child][residue_loc]
        parent_ident = sequence_dict[parent][residue_loc]
        child_pos = residue_dict.get(child_ident)
        parent_pos = residue_dict.get(parent_ident)
        
        if (child_pos != None and parent_pos != None) and (child_pos != parent_pos):
            residue_mut[parent_pos, child_pos] += 1
        elif child_pos is None or parent_pos is None:
            # print(f"No match has found at site {residue_loc + 1} between the child {edge[0]}:{child_ident} and parent {edge[1]}:{parent_ident}")
            continue
         
    return residue_mut 


def WorkerCountMutation(task_queue, result_queue):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        while True:
            try:
                args = task_queue.get(True, 1)
                if args is None:
                    break
                mut_matrix = CountMutations(*args)
                res_loc = args[-1]
                result_queue.put({res_loc: mut_matrix})
            except queue.Empty:
                continue
            except Exception as e:
                # traceback.print_exc()
                result_queue.put((None, f"ERROR0: {e}"))

    
# def residue_table(residues: List[str]) -> Dict[str, int]:
#     residue_pos = [*range(len(residues))]
#     residue_tuples = [*zip(residues, residue_pos)]
#     residue_dict = {k: v for k, v in residue_tuples}
#     residue_dict_flip = {v: k for k, v in residue_tuples}
#     return residue_dict, residue_dict_flip


# def WriteMutMatrix(mut_matrices: Dict[str, str], 
#                    residue_dict_flip: Dict[str, str], 
#                    wd_current: str) -> None: 

#     fn_path = os.path.join(wd_current, "mutation_count_matrices.txt")
#     colname = ["Site", "Parent:Child", "MutCount", "RowIndex", "ColIndex"]
#     accum_mutation_matrix = np.zeros((20, 20), dtype=np.int64)
#     with open(fn_path, "w") as writer:
#         writer.write("\t".join(colname) + "\n")
#         for pos, mut in mut_matrices:
#             accum_mutation_matrix += mut
#             coo = spy.sparse.coo_matrix(mut)
#             data = coo.data 
#             if len(data) == 0:
#                 line = "\t".join((str(pos+1), "-", "-", "-", "-"))
                
#             else:
#                 data_str = ",".join([*map(str, data)])
#                 parent_child = []
#                 for row, col in zip(coo.row, coo.col):
#                     parent = residue_dict_flip[row]
#                     child = residue_dict_flip[col]
#                     parent_child.append(":".join((parent, child)))
#                 parent_child_str = ",".join(parent_child)
#                 row_str = ",".join([*map(str, coo.row)])
#                 col_str = ",".join([*map(str, coo.col)])
#                 line = "\t".join((str(pos+1), parent_child_str, data_str, row_str, col_str))
#             line += "\n"
#             writer.write(line)

#     accum_mut_fnpath = os.path.join(wd_current, "accum_mutation_matrix.txt")
#     with open(accum_mut_fnpath, "w") as writer:
#         col_index = "\t".join(residues)
#         col_index = "\t".join((" ", col_index)) + "\n"
#         writer.write(col_index)
#         for i in range(accum_mutation_matrix.shape[0]):
#             row_str = "\t".join([*map(str, accum_mutation_matrix[i])])
#             row_str = "\t".join((residues[i], row_str))
#             row_str += "\n"
#             writer.write(row_str)

# def WorkWriteMutMatrix(task_queue, result_queue, file_lock, residue_dict_flip, wd_current):

#     fn_path = os.path.join(wd_current, "mutation_count_matrices.txt")
#     colname = ["Site", "Parent:Child", "MutCount", "RowIndex", "ColIndex"]

#     with open(fn_path, "w") as writer:
#         writer.write("\t".join(colname) + "\n")
#         while True:
#             try:
#                 pos, mut = task_queue.get(timeout=1)
#                 coo = spy.sparse.coo_matrix(mut)
#                 data = coo.data 
#                 if len(data) == 0:
#                     line = "\t".join((str(pos+1), "-", "-", "-", "-")) + "\n"
#                 else:
#                     data_str = ",".join(map(str, data))
#                     parent_child = []
#                     for row, col in zip(coo.row, coo.col):
#                         parent = residue_dict_flip[row]
#                         child = residue_dict_flip[col]
#                         parent_child.append(":".join((parent, child)))
#                     parent_child_str = ",".join(parent_child)
#                     row_str = ",".join(map(str, coo.row))
#                     col_str = ",".join(map(str, coo.col))
#                     line = "\t".join((str(pos+1), parent_child_str, data_str, row_str, col_str)) + "\n"
                
#                 # Acquire lock before writing to the file
#                 with file_lock:
#                     writer.write(line)

#                 task_queue.task_done()
#             except queue.Empty:
#                 break
#             except Exception as e:
#                 result_queue.put(f"ERROR at pos {pos}: {e}")


def count_mutations(treefile: str, 
                    sequence_dict: Dict[str, str],
                    residue_dict: Dict[str, str],
                    protein_len: int,
                    options: process_args.Options,
                    ):

    #compile a list of all the child --> parent relationships in the tree
    child_parent_set, root = get_child_parent_relationships(treefile, options.outgroups, options.species_of_interest)
    task_queue = mp.Queue()
    for res_loc in range(protein_len):
        task_queue.put((child_parent_set, sequence_dict, residue_dict, res_loc))

    mut_results = parallel_task_manager.ParallelComputer(WorkerCountMutation, 
                                                            task_queue, 
                                                            options.nProcessAlg)

    mut_results_dict = {}
    for d in mut_results:
        if isinstance(d, dict):
            for k, v in d.items():
                mut_results_dict[k] = v

    return root.strip(), mut_results_dict


