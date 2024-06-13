#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import time
import shutil
import datetime
import scipy as spy 
import numpy as np
from typing import Dict, List, Tuple, Optional, Union
from . import util
from .. import __version__, residues


class FileLogger(object):

    @staticmethod
    def LogFailAndExit(text: str = "", rd1: Optional[str] = None) -> None:
        if text != "": print(text)
        FileLogger.WriteToLog("\nERROR: An error occurred\n" + text, rd1)
        util.Fail()

    @staticmethod         
    def WriteToLog(text: str, rd1: Optional[str] = None, qWithTime: bool = False) -> None:
        prepend = ""
        if qWithTime:
            prepend = str(datetime.datetime.now()).rsplit(".", 1)[0] + " : "
        with open(rd1 + "Log.txt", 'a') as outfile:
            outfile.write(prepend + text)

    @staticmethod
    def StartLog(wd_base0: str, rd1: str, tree_program: Optional[str] = None, sequence_type: Optional[str] = None) -> None:

        FileLogger.WriteToLog( "Started OrthoFinder version " + __version__ + "\n", rd1, True)
        text = "Command Line: " + " ".join(sys.argv) + "\n\n"

        if tree_program:
            text += f"Tree program: {tree_program}\n"

        if sequence_type:
            text += f"Codon model: {sequence_type}\n"

        text += "\nWorkingDirectory_Base: %s\n" % wd_base0#self.wd_base[0]
        FileLogger.WriteToLog(text, rd1)
    

class FileHandler(FileLogger):  

    def __init__(self):
        # self.base_matrix_format = "MUT%07d"
        self.wd_base = []               # Base: blast, species & sequence IDs, species fasta files - should not request this and then write here
        self.wd_current = None          # Location to write out any new files
        self.rd1 = None 
        self.species_ids_corrected = None
        self.gene_of_interest = None
        # to be modified as appropriate
     
    """ ========================================================================================== """
    
    
    def CreateOutputDirectories(self, options, base_dir):

        if options.name == None:
            self.rd1 = util.CreateNewWorkingDirectory(base_dir + "Results_", 
                                                      sequence_type=options.sequence_type,
                                                      extended_filename=options.extended_filename)
        else:
            self.rd1 = util.CreateNewWorkingDirectory(base_dir + "Results_" + options.name, 
                                                      qDate=False,
                                                      sequence_type=options.sequence_type,
                                                      extended_filename=options.extended_filename)
        self.wd_current = self.rd1 + "Real_Phylogeny/"
        os.mkdir(self.wd_current)
        self.wd_base = [self.wd_current]  
        self.StartLog(self.wd_base[0], self.rd1, tree_program=options.tree_program, sequence_type=options.sequence_type)


    def CheckFileCorrespondance(self, gene: str, state_file: str, tree_file:str):

        gene_statfile = os.path.basename(state_file).split(".", 1)[0]
        gene_treefile = os.path.basename(tree_file).split(".", 1)[0]
        if gene != gene_statfile and gene !=  gene_treefile:
            print('ERROR: Gene files have not been correctly processed')
            util.Fail()
    

    def ProcessesNewAln(self, alnDir):
        """
        Process fasta files and return a Directory object with all paths completed.
        """

        alnExtensions = {"aln", "fasta"} # need to add more potential file extensions

        if not os.path.exists(alnDir):
            print("\nDirectory does not exist: %s" % alnDir)
            util.Fail()
        files_in_directory = sorted([f for f in os.listdir(alnDir) if os.path.isfile(os.path.join(alnDir, f))])
        originalALNFilenames = {}
        excludedFiles = []
        for f in files_in_directory:
            if len(f.rsplit(".", 1)) == 2 and f.rsplit(".", 1)[1].lower() in alnExtensions and not f.startswith("._"):
                if len(f.split(".", 1)) == 2:
                    self.gene_of_interest = f.split(".", 1)[0]
                    originalALNFilenames[self.gene_of_interest] = os.path.join(alnDir, f)
            else:
                excludedFiles.append(os.path.join(alnDir, f))

        if len(excludedFiles) != 0:
            print("\nWARNING: Files have been ignored as they don't appear to be ALN files:")
            for f in excludedFiles:
                print(f)
            print("RecFinder expects ALN files to have one of the following extensions: %s" % (", ".join(alnExtensions)), end="\n"*2)
        
        if len(originalALNFilenames) == 0:
            print("\nNo fasta files found in supplied directory: %s" % alnDir)
            util.Fail()

        return originalALNFilenames

    def GetStateFileFN(self):
        if self.wd_current == None: 
            raise Exception("No wd_current")
        statefile = [os.path.join(self.wd_current, file) for file in os.listdir(self.wd_current) if file.endswith(".state")][0]
        return statefile 
    
    def GetTreeFileFN(self):
        if self.wd_current == None: 
            raise Exception("No wd_current")
        treefile = [os.path.join(self.wd_current, file) for file in os.listdir(self.wd_current) if file.endswith(".treefile")][0]
        return treefile 

    def GetIQTreeFileFN(self):
        if self.wd_current == None: 
            raise Exception("No wd_current")
        iqtreefile = [os.path.join(self.wd_current, file) for file in os.listdir(self.wd_current) if file.endswith(".iqtree")][0]
        return iqtreefile 

    def GetMutCountMatricesFN(self):
        if self.GetMutMatrixDir() == None: 
            raise Exception("No Mutation_Matrices directory")
        mut_matrix_file = os.path.join(self.GetMutMatrixDir(), "mutation_matrices.tsv")
        return mut_matrix_file 

    def GetAccumMutCountMatricesFN(self):
        if self.GetMutMatrixDir() == None: 
            raise Exception("No Mutation_Matrices directory")
        mut_matrix_file = os.path.join(self.GetMutMatrixDir(), "accum_mutation_matrix.tsv")
        return mut_matrix_file 

    def GetCombinedDNASeqsFN(self):
        if self.GetInferedSeqsDir() == None: 
            raise Exception("No Infered_Sequences directory")
        combined_seqs_file = os.path.join(self.GetInferedSeqsDir(), f"{self.gene_of_interest}.combined_dna_sequences.aln")
        return combined_seqs_file 

    
    def GetNodeDNASeqsFN(self):
        if self.GetInferedSeqsDir() == None: 
            raise Exception("No Infered_Sequences directory")
        combined_seqs_file = os.path.join(self.GetInferedSeqsDir(), f"{self.gene_of_interest}.node_dna_sequences.aln")
        return combined_seqs_file 

    def GetCombinedProtSeqsFN(self):
        if self.GetInferedSeqsDir() == None: 
            raise Exception("No Infered_Sequences directory")
        combined_seqs_file = os.path.join(self.GetInferedSeqsDir(), f"{self.gene_of_interest}.combined_protein_sequences.aln")
        return combined_seqs_file 

    
    def GetNodeProtSeqsFN(self):
        if self.GetInferedSeqsDir() == None: 
            raise Exception("No Infered_Sequences directory")
        combined_seqs_file = os.path.join(self.GetInferedSeqsDir(), f"{self.gene_of_interest}.node_protein_sequences.aln")
        return combined_seqs_file 

    def GetRecurrenceListFN(self):
        if self.GetRecurrenceStatsDir() == None: 
            raise Exception("No recurrence_Analysis directory")
        recurrence_list_file = os.path.join(self.GetRecurrenceStatsDir(), f"{self.gene_of_interest}.recurrence_list.tsv")
        return recurrence_list_file

    def GetMutMatrixDir(self):
        d = self.rd1 + "Mutation_Matrices/"
        if not os.path.exists(d): os.mkdir(d)
        return d


    def GetInferedSeqsDir(self):
        d = self.rd1 + "Infered_Sequences/"
        if not os.path.exists(d): os.mkdir(d)
        return d

    def GetMCSimulationDir(self):
        d = self.rd1 + "Monte_Carlo_Similation/"
        if not os.path.exists(d): os.mkdir(d)
        return d

    def GetMCSimulationPhyDir(self):
        d = self.GetMCSimulationDir() + "Phy_Files/"
        if not os.path.exists(d): os.mkdir(d)
        return d

    def GetMCSimulationAlnDir(self):
        d = self.GetMCSimulationDir() + "Alignment_Files/"
        if not os.path.exists(d): os.mkdir(d)
        return d

    def GetMCSimulationTreeDir(self):
        d = self.GetMCSimulationDir() + "Tree_File/"
        if not os.path.exists(d): os.mkdir(d)
        return d

    def GetRecurrenceStatsDir(self):
        d = self.rd1 + "Recurrence_Analysis/"
        if not os.path.exists(d): os.mkdir(d)
        return d

    def GetWorkingDirectory1_Read(self):
        if len(self.wd_base) == 0: raise Exception("No wd1")
        return self.wd_base 
        
    def GetWorkingDirectory_Write(self):
        if self.wd_current == None: raise Exception("No wd_current")
        return self.wd_current 
        
    def GetResultsDirectory1(self):
        if self.rd1 == None: raise Exception("No rd1")
        return self.rd1
        

class FileReader(object):

    @staticmethod    
    def ReadAlignment(fn: str) -> tuple[dict: [str, str], int]:
        msa = dict()
        accession = None
        length = None
        seq = ""
        with open(fn, 'r') as infile:
            for line in infile:
                line = line.rstrip()
                if line.startswith(">"):
                    if accession != None:
                        if length != None and len(seq) != length:
                            text  = "ERROR: Sequence length mismatch in MSA: %s & %d" % (length, len(seq))
                            FileLogger.LogFailAndExit(text)
                        msa[accession] = seq
                    accession = line[1:].replace('_', ' ')
                    seq = ""
                else:
                    seq += line
            if accession != None:
                if length != None and len(seq) != length:
                    text = "ERROR: Sequence length mismatch in MSA: %s & %d" % (length, len(seq))
                    FileLogger.LogFailAndExit(text)
                msa[accession] = seq

        return msa, len(seq)
    
    @staticmethod
    def ReadStateFile(fn: str) -> Dict[str, str]:

        node_seq_dict = {}
        with open(fn, 'r') as reader:
            for line in reader:
                if "#" in line or "Site" in line:
                    continue
                line = line.replace("\n", "").strip().split("\t", 3)
                node, state = line[0], line[2]
                if node not in node_seq_dict:
                    node_seq_dict[node] = state
                elif node in node_seq_dict:
                    node_seq_dict[node] += state

        return node_seq_dict
    @staticmethod
    def ReadMutMatrixFile(mut_matrix_path: str):

        matrix_dict ={}
        mut_matrix = np.zeros((20, 20))
        with open(mut_matrix_path, 'r') as reader:
            for i, line in enumerate(reader):
                if "Site" in line:
                    continue
                line = line.replace("\n", "").strip().split()
                data = line[2]
                if data == "-":
                    matrix_dict[i] = mut_matrix

                else:
                    data = np.array([*map(np.int64, data.split(","))])
                    row = np.array([*map(np.int64, line[3].split(","))])
                    col = np.array([*map(np.int64, line[4].split(","))])#
                    coo = spy.sparse.coo_matrix((data, (row, col)), shape=(20, 20))
                    matrix = coo.toarray()
                    matrix_dict[i] = matrix

        return matrix_dict

    @staticmethod
    def ReadIQTreeFile(iqtree_path: str):
        with open(iqtree_path) as reader:
            for line in reader:
                if "Best-fit model" in line:
                    best_evolution_model = line.replace("\n", "").strip().split(": ")[-1]
        return best_evolution_model


    @staticmethod
    def ReadPhyFile(phy_path: str):

        phy_dict = {}
        with open(phy_path) as reader:
            for line in reader:
                line = line.replace("\n", "").strip().split()
                try: 
                    val0 = int(line[0])
                    val1 = int(line[1])
                except:
                    phy_dict[line[0]] = line[1]

        return phy_dict


class FileWriter(object):
    
    @staticmethod
    def WriteMutMatrix(mut_matrices: Dict[int, np.ndarray], 
                        residue_dict_flip: Dict[str, str], 
                        mut_matrix_dir: str) -> None: 

        fn_path = os.path.join(mut_matrix_dir, "mutation_matrices.tsv")
        colname = ["Site", "Parent:Child", "MutCount", "RowIndex", "ColIndex"]
        accum_mutation_matrix = np.zeros((20, 20), dtype=np.int64)
        with open(fn_path, "w") as writer:
            writer.write("\t".join(colname) + "\n")
            for pos, mut in mut_matrices.items():
                accum_mutation_matrix += mut
                coo = spy.sparse.coo_matrix(mut)
                data = coo.data 
                if len(data) == 0:
                    line = "\t".join((str(pos+1), "-", "-", "-", "-"))
                    
                else:
                    data_str = ",".join(map(str, data))
                    parent_child = []
                    for row, col in zip(coo.row, coo.col):
                        parent = residue_dict_flip[row]
                        child = residue_dict_flip[col]
                        parent_child.append(":".join((parent, child)))
                    parent_child_str = ",".join(parent_child)
                    row_str = ",".join(map(str, coo.row))
                    col_str = ",".join(map(str, coo.col))
                    line = "\t".join((str(pos+1), parent_child_str, data_str, row_str, col_str))
                line += "\n"
                writer.write(line)

        accum_mut_fnpath = os.path.join(mut_matrix_dir, "accum_mutation_matrix.tsv")
        with open(accum_mut_fnpath, "w") as writer:
            col_index = "\t".join(residues)
            col_index = "\t".join((" ", col_index)) + "\n"
            writer.write(col_index)
            for i in range(accum_mutation_matrix.shape[0]):
                row_str = "\t".join(map(str, accum_mutation_matrix[i]))
                row_str = "\t".join((residues[i], row_str))
                row_str += "\n"
                writer.write(row_str)


    @staticmethod
    def WriteSeqsToAln(seqs_dict: Dict[str, str], outFilename: str) -> None:
        with open(outFilename, 'w') as writer:
            for seq in seqs_dict:
                writer.write(">%s\n" % seq)
                writer.write(seqs_dict[seq] + "\n")
                writer.write("\n")

    @staticmethod
    def WriteRecurrenceList(recurrence_list: List[List[Union[str, int, float]]],
                            outFilename: str) -> None:

        colname = ['Site', 'Parent', 'Child', 'Recurrence', 
                   "Reversion", "P-Value", "OtherSiteSubs",  "SiteComposition"]
        with open(outFilename, "w") as writer:
            writer.write("\t".join(colname) + "\n")
            for rec_list in recurrence_list:
                rec_list[0] += 1
                writer.write("\t".join(map(str, rec_list)) + "\n")



        