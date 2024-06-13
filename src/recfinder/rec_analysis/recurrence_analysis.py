#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 09:52:22 2023

This script will get a dataframe of recurrence of repeated mutations
@author: Lizzie
"""
import pandas as pd
from pathlib import Path
from ..utils import files, util, process_args, parallel_task_manager
from typing import Optional, Dict, List, Tuple, Union
from .. import residues
import os
import numpy as np
import ast
import dendropy
from ..utils import process_args
import scipy as spy
import multiprocessing as mp
import warnings
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
try:
    import queue
except ImportError:
    import Queue as queue     
import cProfile


def profile(func):
    def wrapper(*args, **kwargs):
        profile = cProfile.Profile()
        profile.enable()
        result = func(*args, **kwargs)
        profile.disable()
        profile.print_stats()
        return result
    return wrapper


def GetrecurrenceList(rec_loc: int,
                      mut_matrix: np.ndarray, 
                      ident: Tuple[str, ...]) -> List[Union[int, float, str]]:

    recurrence_list = []
    rows, cols = np.where(mut_matrix > 1) # why not m > 0?
    for i in range(0, len(rows)):
        r = rows[i]
        c = cols[i]
        recurrence = mut_matrix[r, c]
        flipflop = mut_matrix[c, r]
        parent = residues[r]
        child = residues[c]

        parent_percent = (ident.count(parent) / len(ident)) * 100
        child_percent = (ident.count(child) / len(ident)) * 100
        
        recurrence_list.append([rec_loc, parent, child, recurrence, flipflop])
    
    return recurrence_list

def WorkerrecurrenceList(task_queue: mp.Queue, result_queue: mp.Queue):

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        while True:
            try:
                args = task_queue.get(True, 1)
                if args is None:
                    break
                recurrence_list = GetrecurrenceList(*args)
                rec_loc = args[0]
                result_queue.put((rec_loc, recurrence_list))
            except queue.Empty:
                continue
            except Exception as e:
                # traceback.print_exc()
                result_queue.put((None, f"ERROR: {e}"))


def recurrence_list(mut_matrices:  Dict[int, np.ndarray], 
                    extant_seq: Dict[str, str], 
                    outgroups: List[str],
                    nthreads: int) -> List[List[Union[str, int]]]:

    sequence_dict = {species: seq for species, seq in extant_seq.items() if species not in outgroups}
    ident_dict = {}
    for rec_loc, res in enumerate(zip(*sequence_dict.values())):
        ident_dict[rec_loc] = res
    
    task_queue = mp.Queue()
    for rec_loc, mut_matrix in mut_matrices.items():
        task_queue.put((rec_loc, mut_matrix, ident_dict[rec_loc]))

    recurrence_results = parallel_task_manager.ParallelComputer(WorkerrecurrenceList, 
                                                                task_queue, 
                                                                nthreads)
    recurrence_list = []
    for rec_loc, rec_list in recurrence_results:
        if len(rec_list) == 0:
            continue 
        recurrence_list.extend(rec_list)

    return recurrence_list, ident_dict


def WorkerComputePValues(task_queue, result_queue): 
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        while True:
            try:
                args = task_queue.get(True, 1)
                if args is None:
                    break
                count_greater, p_value = ComputePValues(*args)
                result_queue.put((*args, count_greater, p_value))
            except queue.Empty:
                continue
            except Exception as e:
                # traceback.print_exc()
                result_queue.put((None, f"ERROR: {e}"))

        
def compute_p_values(mcs_results: List[Tuple[str, str, Dict[int, np.ndarray]]],
                     recurrence_list: List[List[Union[str, int, float]]],
                     residue_dict: Dict[str, int],
                     nthreads: int):


    def compute_p_value_for_rec_list(rec_list):
        rec_loc, parent, child, recurrence = rec_list[:4]
        r = residue_dict[parent]
        c = residue_dict[child]

        count_greater = 0

        for mcs_result in mcs_results:
            if len(mcs_result) < 3:
                continue

            if all(np.all(array == 0) for array in mcs_result[-1].values()):
                continue

            m = mcs_result[-1].get(rec_loc)

            if m is None:
                continue

            mcs_rec = m[r, c]
            if mcs_rec > recurrence:
                count_greater += 1

        return count_greater / 1000

    recurrence_list_pvalue = []
    with ThreadPoolExecutor(max_workers=nthreads) as executor:
        future_to_rec_list = {executor.submit(compute_p_value_for_rec_list, rec_list): rec_list for rec_list in recurrence_list}
        
        for future in as_completed(future_to_rec_list):
            rec_list = future_to_rec_list[future]
            try:
                p_value = future.result()
                rec_list.append(p_value)
            except Exception as exc:
                traceback.print_exc()
                print(f"Generated an exception: {exc}")
                rec_list.append(None)

            print(rec_list)
            recurrence_list_pvalue.append(rec_list)

    recurrence_list_pvalue = sorted(recurrence_list_pvalue, key=lambda x: (x[-1], -x[-3]))
    return recurrence_list_pvalue
    
    
    
    
    
    # mcs_mut_matrices = []
    # for mcs_result in mcs_results:
    #     if len(mcs_result) < 3:
    #         continue
    #     check_all_zero = all(np.all(array == 0) for array in mcs_result[-1].values())
    #     if check_all_zero:
    #         continue 
    #     mcs_mut_matrices.append(mcs_result[-1])

    # task_queue = mp.Queue()
    # for rec_list in recurrence_list:
    #     rec_loc, parent, child, recurrence = rec_list[:4]
    #     r = residue_dict[parent]
    #     c = residue_dict[child]

    #     count_greater = 0
    #     for mcs_result in mcs_results:
    #         if len(mcs_result) < 3:
    #             continue
    #         check_all_zero = all(np.all(array == 0) for array in mcs_result[-1].values())
            
    #         if check_all_zero:
    #             continue 
    #         m = mcs_result[rec_loc]
    #         mcs_rec = m[r, c]
    #         if mcs_rec > recurrence:
    #             count_greater += 1

    #     p_value = count_greater / 1000
        # mut_matrix = mcs_results[-1][rec_loc]
        # task_queue.put((rec_loc, r, c, recurrence, mcs_mut_matrices))

    # for mcs_result in mcs_results:
    #     if len(mcs_result) > 2:
    #         mcs_mut_matrices = mc_result[-1] 
    #         task_queue.put(mcs_mut_matrices)

    # pvalue_results = parallel_task_manager.ParallelComputer(WorkerComputePValues, 
    #                                                             task_queue, 
    #                                                             nthreads)

    # for item in pvalue_results:
    #     print(item)

#     for site in dict_keys_tuple:
#         x = list(site)
#         k = "".join(x)
#         count_dict[k] = []

#     for index, residue_path in enumerate(matrix_paths):
#         print(index)
#         d, sequence_lines = get_mutation_matrix(residue_path)
#         #print(d[1].shape)
#         for index, row in gene_df.iterrows():
#             parent = row['Parent']
#             r = residues.index(parent)
#             child = row['Child']
#             c = residues.index(child)

#             res = row['Site']
#             mut_matrix = d[res]
#             recurrence = mut_matrix[r,c]

#             identifier = str(res) + parent + child
#             count_dict[identifier].append(recurrence)

#     output_list = []
#     for identifier, count_list in count_dict.items():
#         res = int(identifier[:-2])
#         res_idents = identifier[len(identifier)-2:]
#         parent = res_idents[0]
#         child = res_idents[1]
#         row_of_interest = gene_df[(gene_df['Site'] == res) & (gene_df['Parent'] == parent) & (gene_df['Child'] == child)]
#         count_to_beat = row_of_interest['recurrence'].item()

#         count_greater = len([i for i in count_list if i >= count_to_beat])
#         p_value = count_greater/1000
#         output_row = [gene, res, parent, child, count_to_beat, count_greater, p_value]
#         output_list.append(output_row)

#     out_df = pd.DataFrame(output_list, columns = ['Gene', 'Site', 'Parent', 'Child', 'recurrence', 'Count_>=_recurrence', 'p-value'])
#     output_df_path = script_directory + '/' + gene + '_p_values.csv'
#     out_df.to_csv(output_df_path, index = False)

#     output_dict_path = script_directory + '/' + gene + '_count_dict.json'

#     with open(output_dict_path, 'w') as f:
#         for res, counts in count_dict.items():
#             f.write(">" + res + '\n' + str(counts) + '\n' + '\n')

#     paths_analysed_path = script_directory + '/' + gene + '_paths_analysed.txt'
#     with open(paths_analysed_path, 'w') as f:
#         for p in matrix_paths:
#             f.write(p + '\n')