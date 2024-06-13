#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 10 15:19:57 2022

@author: Lizzie (elizabeth.robbins@wlofson.ox.ac.uk)

This script counts amino acid mutations that have occured throughout a phylogenetic tree taking into account a ancestral state reconstruction (from IQTREE).

Notes:
    State files only contain the ancestral state reconstruction sequences (not extant sequences).
    Therefore, in order to get all child --> parent residue identities the sequence alignments also need to be read in.

"""

#IMPORT LIBRARIES

from __future__ import absolute_import

# first import parallel task manager to minimise RAM overhead for small processes
import multiprocessing as mp                   
import platform                                
import sys                                     

if __name__ == "__main__":
    if platform.system() == "Darwin":
        mp.set_start_method('fork')


import os
import csv
import time
import cProfile
import dendropy
import pandas as pd
import numpy as np
from tqdm import tqdm
import json
from typing import Optional, Dict, List, Tuple, Union
import multiprocessing as mp
from collections import Counter
import scipy as spy
import warnings
import traceback
from .. import residues
from .. import  __version__
from . import run_commands
from ..utils import files, util, parallel_task_manager, process_args, timeit
from ..rec_analysis import count_aa_mutations, recurrence_analysis
try:
    import queue
except ImportError:
    import Queue as queue 


max_int = sys.maxsize
ok = False
while not ok:
    try:
        csv.field_size_limit(max_int)
        ok = True
    except OverflowError:
        max_int = int(max_int/10)
sys.setrecursionlimit(10**6)


def residue_table() -> Dict[str, int]:

    residue_pos = [*range(len(residues))]
    residue_tuples = [*zip(residues, residue_pos)]
    residue_dict = {k: v for k, v in residue_tuples}
    residue_dict_flip = {v: k for k, v in residue_tuples}

    return residue_dict, residue_dict_flip

def GetGeneTreeBuildCommands(alignment_file_list: list[str], 
                              wd_current: str, 
                              nthreads: int,
                              evolution_model: str,
                              output_prefix: Optional[str] = None,
                              sequence_type: Optional[str] = None,
                              gene_tree: Optional[str] = None,
                              partial_constraints: Optional[str] = None) -> list[str]:
    commands = []
    for alignment_file in alignment_file_list:
        identifier = os.path.basename(alignment_file).rsplit(".", 1)[0]
        if not output_prefix:
            pre = wd_current + identifier
        else:
            pre = output_prefix + identifier
        command = run_commands.GetGeneTreeBuildCommand(alignment_file, 
                                                nthreads,
                                                evolution_model,
                                                pre, 
                                                sequence_type=sequence_type,
                                                gene_tree=gene_tree, 
                                                partial_constraints=partial_constraints)
        commands.append(command)
    return commands


def WorkerReadPhyFile(task_queue, result_queue):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        while True:
            try:
                phy_path = task_queue.get(timeout=1)
                if phy_path is None:  # Sentinel to stop the worker
                    break
                try:
                    result = files.FileReader.ReadPhyFile(phy_path)
                    result_queue.put((phy_path, result))
                except Exception as e:
                    # traceback.print_exc()
                    result_queue.put((phy_path, f"ERROR: {e}"))
            except queue.Empty:
                continue  # Keep checking for tasks or sentinel

def WorkerCountAAMutation(task_queue, result_queue):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        while True:
            try:
                args = task_queue.get(True, 1)
                if args is None:
                    break
                root_node, mut_matrices = count_aa_mutations.count_mutations(*args)
                fn = os.path.basename(args[0]).split(".")[0]
                
                # mut_matrices_3d = np.stack([m for _, m in mut_matrices])
                # print(mut_matrices_3d.shape)
                result_queue.put((fn, root_node, mut_matrices))
            except queue.Empty:
                continue
            except Exception as e:
                # traceback.print_exc()
                result_queue.put((None, f"ERROR1: {e}"))

def combine_mut_matrix_and_recurrence_list(mut_matrices: Dict[int, np.ndarray],
                                           recurrence_list: List[List[Union[str, int]]],
                                           ident_dict: Dict[int, List[str]],
                                           residue_dict_flip: Dict[int, str]) -> List[List[Union[str, int]]]:

    for rec_list in recurrence_list:
        res_loc = rec_list[0]
        mut = mut_matrices[res_loc]
        coo = spy.sparse.coo_matrix(mut)
        data, row, col = coo.data, coo.row, coo.col
        sorted_indices = np.argsort(data)
        data = data[sorted_indices]
        row = row[sorted_indices]
        col = col[sorted_indices]

        data_str_list = [*(map(str, data))]
        parent_child = []
        for row, col in zip(row, col):
            parent = residue_dict_flip[row]
            child = residue_dict_flip[col]
            parent_child.append(">".join((parent, child)))
        parent_child_data = [*zip(parent_child, data_str_list)]
        parent_child_data_str = ",".join([":".join(pcd) for pcd in parent_child_data])
        rec_list.append(parent_child_data_str)

        res_freq = [*Counter(ident_dict[res_loc]).items()]
        res_freq = sorted(res_freq, reverse=True, key=lambda x: x[1])
        res_freq_str = ",".join([":".join((res, str(freq))) for res, freq in res_freq])
        rec_list.append(res_freq_str)

    return recurrence_list


@timeit.timeit
def main(args: Optional[list[str]] = None):
    
    try:
        if not args:
            args = sys.argv[1:]
        
        ## ------------------------------------------------------------------------------------
        # create a pool of processes as early in execution as possible so that the memory footprint is low
        ptm_initialised = parallel_task_manager.ParallelTaskManager_singleton()
        options, alnDir, resultsDir_nonDefault = process_args.ProcessArgs(args)
        
        # create output directory
        base_dir = resultsDir_nonDefault if resultsDir_nonDefault != None else alnDir + "RecFinder/"
        if not os.path.exists(base_dir): 
            os.mkdir(base_dir)
        
        # Initiate the file instances 
        filehandler = files.FileHandler()
        filereader = files.FileReader()
        filewriter = files.FileWriter()
        
        # Create the output directory
        filehandler.CreateOutputDirectories(options, base_dir)

        # Extract the new alignment file information
        aln_path_dict = filehandler.ProcessesNewAln(alnDir)
        gene = [*aln_path_dict.keys()][0]
        print(f"Analysing: {gene}")
        results_dir = filehandler.GetResultsDirectory1()
        print(f"Results Directory: {results_dir}", end="\n"*2)

        ## ------------------------------------------------------------------------------------------
        step1_info = "Step1: Analysing Real Phylogeny"
        print(step1_info)
        print("="*len(step1_info))
        print(f"Results Directory: {filehandler.wd_current}", end="\n"*2)

        if not options.usr_state or not options.usr_tree:

            print_info1 = "Ran IQ-TREE to build the gene trees"
            commands = GetGeneTreeBuildCommands([*aln_path_dict.values()], 
                                                filehandler.wd_current, 
                                                options.nMSATree,
                                                evolution_model = options.evolution_model,
                                                sequence_type = options.sequence_type,
                                                gene_tree = options.gene_tree,
                                                partial_constraints = options.partial_constraints)
            run_commands.RunCommand(options.nMSATree, commands, print_info1, 
                                    filehandler.wd_current, delete_files=True,
                                    files_to_keep=["state", "treefile", "iqtree"])          
            util.PrintTime("Done gene tree building\n")
            
            statefile = filehandler.GetStateFileFN() 
            treefile = filehandler.GetTreeFileFN() 

        else:
            statefile = options.usr_state
            treefile = options.usr_tree

        # statefile = r"./ExampleData/RecFinder/Results_Jun11_1/Real_Phylogeny/28_species_rbcL.nuc.state" # r'./rbcL.nuc.aln.state'
        # treefile = r"./ExampleData/RecFinder/Results_Jun11_1/Real_Phylogeny/28_species_rbcL.nuc.treefile" # r'./rbcL.nuc.aln.treefile'
        
        # Defines species in the tree of interest
        outgroup_mrca = options.outgroups
        mrca_of_interest = options.species_of_interest

        filehandler.CheckFileCorrespondance(gene, statefile, treefile)

        statefile_dict = {gene: statefile}
        treefile_dict = {gene: treefile}
        alignment_file = aln_path_dict[gene]

        alignment_dict, alignment_len = filereader.ReadAlignment(alignment_file)
        print(f"Number of extant species found: {len(alignment_dict)}")
        print(f"Length of the {gene} alignment: {alignment_len}")
            
        # codon_table = import_codon(options.sequence_type)     #NCBI's genetic code 11 (for plant plastids)
        residue_dict, residue_dict_flip = residue_table()
        
        # get the ancestral state reconstruction sequences for each internal node
        node_seq_dict = filereader.ReadStateFile(statefile)
        
        # merge alignment dictionary and state dictionary (have all sequences in one place)
        combined_seq_dict = {k: v for d in (node_seq_dict, alignment_dict) for k, v in d.items()}
        
        if options.sequence_type == "AA": # No translation is needed when the input sequence type is amino acid
            node_prot_seqs_fn = filehandler.GetNodeProtSeqsFN()
            combined_prot_seqs_fn = filehandler.GetCombinedProtSeqsFN()

            filewriter.WriteSeqsToAln(node_seq_dict, node_prot_seqs_fn)
            filewriter.WriteSeqsToAln(combined_seq_dict, combined_prot_seqs_fn)
            
            protein_len = len([*node_seq_dict.values()][0])

        else:
            node_dna_seqs_fn = filehandler.GetNodeDNASeqsFN()
            combined_dna_seqs_fn = filehandler.GetCombinedDNASeqsFN()

            filewriter.WriteSeqsToAln(node_seq_dict, node_dna_seqs_fn)
            filewriter.WriteSeqsToAln(combined_seq_dict, combined_dna_seqs_fn)

            combined_prot_seqs_dict, protein_len = util.GetSeqsDict(combined_seq_dict, options.sequence_type)
            node_prot_seqs_dict, _ = util.GetSeqsDict(node_seq_dict, options.sequence_type)
    
            node_prot_seqs_fn = filehandler.GetNodeProtSeqsFN()
            combined_prot_seqs_fn = filehandler.GetCombinedProtSeqsFN()
            
            filewriter.WriteSeqsToAln(node_prot_seqs_dict, node_prot_seqs_fn)
            filewriter.WriteSeqsToAln(combined_prot_seqs_dict, combined_prot_seqs_fn)
        
        # Create mutation matrices at each site and count the substitutions for the real phylogeny
        root_node, mut_matrices = count_aa_mutations.count_mutations(treefile, 
                                                                    combined_prot_seqs_dict,
                                                                    residue_dict, 
                                                                    protein_len, 
                                                                    options)

        print("Root is: ", root_node, end="\n"*2)
        filewriter.WriteMutMatrix(mut_matrices, residue_dict_flip, filehandler.GetMutMatrixDir())
        
        recurrence_list, ident_dict = recurrence_analysis.recurrence_list(mut_matrices, 
                                                                            combined_prot_seqs_dict, 
                                                                            options.outgroups,
                                                                            options.nProcessAlg)
        
        # ## -------------------------------------------------------------------------------------------------
        step2_info = f"\nStep2: Simulating Sequence Evolution with {options.nalign} simulated phylogenies"
        print(step2_info)
        print("="*len(step2_info))
        print(f"Results Directory: {filehandler.GetMCSimulationDir()}", end="\n"*2)
        
        iqtreefile = filehandler.GetIQTreeFileFN()
        # iqtreefile = r"./ExampleData/RecFinder/Results_Jun09/RealPhylogeny/28_species_rbcL.nuc.iqtree"
        best_evolution_model = filereader.ReadIQTreeFile(iqtreefile)
        
        mcs_phyDir = filehandler.GetMCSimulationPhyDir() 
        identifier = "rooted_" + gene + "_alisim"
        output_prefix = os.path.join(mcs_phyDir, identifier)
        if options.sequence_type == "AA":
            root_node = ",".join((combined_prot_seqs_file, root_node))
        else:
            root_node = ",".join((combined_dna_seqs_fn, root_node))
        mc_commands = run_commands.GetMCsimulationCommand(output_prefix, 
                                                        options.nMSATree,
                                                        best_evolution_model, 
                                                        treefile, 
                                                        root_node, 
                                                        options.nalign)
        
        print_info2 = "Ran Monte-Carlo Simulation"
        run_commands.RunCommand(options.nMSATree, mc_commands, print_info2,
                                mcs_phyDir, 
                                delete_files=True,
                                files_to_keep=["phy"])          
        util.PrintTime("Done Monte-Carlo Simulation\n")

        ## -------------------------------------------------------------------------------------------------------------------
        step3_info = f"\nStep3: Analysing the {options.nalign} simulated phylogenies"
        print(step3_info)
        print("="*len(step3_info))
        print(f"Results Directory: {filehandler.GetMCSimulationDir()}", end="\n"*2)

        mcs_alnDir = filehandler.GetMCSimulationAlnDir()
        mcs_phyDir = filehandler.GetMCSimulationPhyDir()  
        mcs_treeDir = filehandler.GetMCSimulationTreeDir()

        print_info3 = "Ran IQ-TREE to build gene trees for the simulated phylogeny"
        mcs_aln_flist = [os.path.join(mcs_phyDir, file) for file in os.listdir(mcs_phyDir)]
        commands = GetGeneTreeBuildCommands(mcs_aln_flist, 
                                            mcs_treeDir, 
                                            options.nMSATree,
                                            output_prefix = mcs_treeDir,
                                            evolution_model = best_evolution_model,
                                            sequence_type = options.sequence_type,
                                            gene_tree = treefile,
                                            partial_constraints = options.partial_constraints)

        run_commands.RunCommand(options.nMSATree, commands, print_info3, mcs_treeDir, 
                                convert_statefiles=True, 
                                sequence_type=options.sequence_type,
                                delete_files=True,
                                files_to_keep=["treefile"])  
            
        util.PrintTime("Done gene tree building\n")
        

        ## ----------------------------------------------------------------------------------------------------------------
        # mcs_alnDir = r"./ExampleData/RecFinder/Results_Jun11_1/Monte_Carlo_Similation/Alignment_Files"
        # mcs_treeDir = r"./ExampleData/RecFinder/Results_Jun11_1/Monte_Carlo_Similation/Tree_File"
        
        recurrenceDir = filehandler.GetRecurrenceListFN()

        step4_info = f"Step4: Creating a recurrence list for the simulated phylogeny"
        print(step4_info)
        print("="*len(step4_info))
        print(f"Results Directory: {recurrenceDir}", end="\n"*2)
        
        task_queue = mp.Queue()
        for file in os.listdir(mcs_alnDir):
            file_path = os.path.join(mcs_alnDir, file)
            mcs_combined_prot_seqs_dict, mcs_protein_len = filereader.ReadAlignment(file_path)
            mcs_treefile = os.path.join(mcs_treeDir, file.split(".")[0] + ".treefile")
            if not os.path.exists(mcs_treefile):
                print("ERROR: Simulated treefile does not exist!")
                util.Fail()
            
            task_queue.put((mcs_treefile, 
                            mcs_combined_prot_seqs_dict, 
                            residue_dict,
                            mcs_protein_len,
                            options))

        # Create mutation matrices at each site and count the substitutions for the simulated phylogeny
        mcs_results = parallel_task_manager.ParallelComputer(WorkerCountAAMutation, 
                                                            task_queue, 
                                                            options.nProcessAlg)


        recurrence_list_pvalue = recurrence_analysis.compute_p_values(mcs_results, recurrence_list, residue_dict, options.nProcessAlg)
        

        recurrence_list_updated = combine_mut_matrix_and_recurrence_list(mut_matrices,
                                                                recurrence_list_pvalue,
                                                                ident_dict,
                                                                residue_dict_flip)

        filewriter.WriteRecurrenceList(recurrence_list_updated, recurrenceDir)
        
        
        d_results = os.path.normpath(filehandler.GetResultsDirectory1()) + os.path.sep
        print("\nResults:\n    %s" % d_results)
        util.PrintCitation(d_results)
        files.FileLogger.WriteToLog("RecFinder run completed\n", d_results, True)

    except Exception as e:
        print(str(e))
        parallel_task_manager.print_traceback(e)
        ptm = parallel_task_manager.ParallelTaskManager_singleton()
        ptm.Stop()
        raise

    except KeyboardInterrupt:
        print("\nProgram terminated by user.")
        sys.exit(1)
    
    finally:
        ptm = parallel_task_manager.ParallelTaskManager_singleton()
        ptm.Stop()

        
if __name__ == "__main__":

    args = sys.argv[1:]
    main(args)

