# -*- coding: utf-8 -*-
#
import gzip
import os
import sys
import time
import copy
import numpy as np
import subprocess
import datetime
try: 
    import queue
except ImportError:
    import Queue as queue
from collections import namedtuple
from ..citation import citation, print_citation
import contextlib
from . import parallel_task_manager, files
from typing import List, Dict, Tuple
from .. import genetic_codes
from importlib import resources as impresources



def ImportCodon(code_name: str):
    
    code_fn = code_name.lower() + ".txt"
    codon_table = {}
    try:
        codon_file = (impresources.files(genetic_codes) / code_fn)

        with codon_file.open("rt") as reader:  # "rt" as text file with universal newlines
            for line in reader:
                codon, letter = line.replace("\n", "").rstrip().split()
                codon_table[codon] = letter
            
    except FileNotFoundError as e:
            print(f"File not found: {e.filename}")
    
    return codon_table


def Translate(seq: str, table: dict) -> str:

    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            if codon not in table:
                protein += '-'
            else:
                protein += table[codon]

    return protein


def GetSeqsDict(dna_seq_dict: Dict[str, str], sequence_type: str):

    codon_table = ImportCodon(sequence_type)     #NCBI's genetic code 11 (for plant plastids)

    #else you need to translate the nucleotide sequence into a protein sequence
    prot_sequence_dict = {}
    for node, seq in dna_seq_dict.items():
        prot_seq = Translate(seq, codon_table)
        prot_sequence_dict[node] = prot_seq
    protein_len = len(prot_seq)

    return prot_sequence_dict, protein_len


def PrintNoNewLine(text):
    parallel_task_manager.PrintNoNewLine(text)

def PrintTime(message):
    parallel_task_manager.PrintTime(message)   
      
"""
Directory and file management
-------------------------------------------------------------------------------
"""               
               
def GetDirectoryName(baseDirName, i, sequence_type, extended_filename):
    
    if not extended_filename:
        if i == 0:
            return baseDirName + os.sep
        else:
            return baseDirName + ("_%d" % i) + os.sep
    else:
        sequence_type = sequence_type if sequence_type else ""
        extension = sequence_type

        if i == 0:
            return baseDirName + "_" + extension + os.sep 
        else:
            return baseDirName + ("_%d" % i) + "_" + extension + os.sep

"""Call GetNameForNewWorkingDirectory before a call to CreateNewWorkingDirectory to find out what directory will be created"""
def CreateNewWorkingDirectory(baseDirectoryName, 
                              qDate=True,
                              sequence_type=None, 
                              extended_filename=False):

    dateStr = datetime.date.today().strftime("%b%d") if qDate else ""
    iAppend = 0
    newDirectoryName = GetDirectoryName(baseDirectoryName + dateStr, 
                                        iAppend,
                                        sequence_type,
                                        extended_filename)
    while os.path.exists(newDirectoryName):
        iAppend += 1
        newDirectoryName = GetDirectoryName(baseDirectoryName + dateStr, 
                                            iAppend, 
                                            sequence_type,
                                            extended_filename)
    os.mkdir(newDirectoryName)
    return newDirectoryName

 
def Success():
    parallel_task_manager.Success()
   
def Fail():
    parallel_task_manager.Fail()
    
"""
IDExtractor
-------------------------------------------------------------------------------
"""

def GetIDPairFromString(line):
    return list(map(int, line.split("_")))

class IDExtractor(object):
    """IDExtractor deals with the fact that for different datasets a user will
    want to extract a unique sequence ID from the fasta file accessions uin different 
    ways."""
    def GetIDToNameDict(self):
        raise NotImplementedError("Should not be implemented")

class FullAccession(IDExtractor):
    def __init__(self, idsFilename):
        # only want the first part and nothing else (easy!)
        self.idToNameDict = dict()
        with open(idsFilename, 'r') as idsFile:
            for line in idsFile:
                line = line.rstrip()
                if not line: continue
#                if line.startswith("#"): continue
                id, accession = line.split(": ", 1)
                id = id.replace("#", "")
                id = id.strip()
                # Replace problematic characters
                accession = accession.replace(":", "_").replace(",", "_").replace("(", "_").replace(")", "_") #.replace(".", "_")
                if id in self.idToNameDict:
                    raise RuntimeError("ERROR: A duplicate id was found in the fasta files: % s" % id)
                self.idToNameDict[id] = accession 
                
    def GetIDToNameDict(self):
        return self.idToNameDict
                
class FirstWordExtractor(IDExtractor):
    def __init__(self, idsFilename):
        # only want the first part and nothing else (easy!)
        self.idToNameDict = dict()
        accs_in_species = []
        with open(idsFilename, 'r') as idsFile:
            for line in idsFile:
                id, rest = line.split(": ", 1)
                accession = rest.split(None, 1)[0]
                iSp = int(id.split("_")[0])
                while len(accs_in_species) < iSp + 1:
                    accs_in_species.append(set())
                # Replace problematic characters
                accession = accession.replace(":", "_").replace(",", "_").replace("(", "_").replace(")", "_") #.replace(".", "_")
                # Only ensure the accessions are unique within the species, there's no need to ensure global uniqueness
                if accession in accs_in_species[iSp]:
                    raise RuntimeError("A duplicate accession was found using just first part: % s" % accession)
                accs_in_species[iSp].add(accession)
                if id in self.idToNameDict:
                    raise RuntimeError("ERROR: A duplicate id was found in the fasta files: % s" % id)
                self.idToNameDict[id] = accession     
                
    def GetIDToNameDict(self):
        return self.idToNameDict

def HaveSupportValues(speciesTreeFN_ids):
    qHaveSupport = False
    try:
        tree.Tree(speciesTreeFN_ids, format=2)
        qHaveSupport = True
    except:
        pass
    return qHaveSupport

def RenameTreeTaxa(treeFN_or_tree, newTreeFilename, idsMap, qSupport, qFixNegatives=False, inFormat=None, label=None, qViaCopy=False):
    """
    qViaCopy - create a copy of the tree and edit this copy. I.e. don't make changes to the original 
    """
    if label != None: qSupport = False
    qHaveSupport = False
    try:
        if type(treeFN_or_tree) is tree.TreeNode:
            if qViaCopy:
                t = treeFN_or_tree.copy("newick")
            else:
                t = treeFN_or_tree
        else:
            qHaveSupport = False
            if inFormat == None:
                try:
                    t = tree.Tree(treeFN_or_tree, format=2)
                    qHaveSupport = True
                except:
                    t = tree.Tree(treeFN_or_tree)
            else:
                t = tree.Tree(treeFN_or_tree, format=inFormat)
        for node in t.get_leaves():
            node.name = idsMap[node.name]
        if qFixNegatives:
            tree_length = sum([n.dist for n in t.traverse() if n != t])
            sliver = tree_length * 1e-6
        iNode = 1
        for n in t.traverse():
            if qFixNegatives and n.dist < 0.0: n.dist = sliver
            if label != None:
                if (not n.is_leaf()) and (not n.is_root()):
                    n.name = label + ("%d" % iNode)
                    iNode += 1
        if label != None:
            with open(newTreeFilename, 'w') as outfile:
                outfile.write(t.write(format=3)[:-1] + label + "0;")  # internal + terminal branch lengths, leaf names, node names. (tree library won't label root node)
        elif t.name == "N0" or t.name == "n0":
            with open(newTreeFilename, 'w') as outfile:
                outfile.write(t.write(format=3)[:-1] + t.name + ";")  # internal + terminal branch lengths, leaf names, node names. (tree library won't label root node)
        else:
            if qSupport or qHaveSupport:
                t.write(outfile = newTreeFilename, format=2)  
            else:
                t.write(outfile = newTreeFilename, format=5)  
    except:
        pass
    
"""
Find results of previous run    
-------------------------------------------------------------------------------
"""

def GetSpeciesDirectory():
    # Confirms all required Sequence files and BLAST etc are present
    pass

def WriteCitation(d):
    with open(d + "Citation.txt", 'w') as outfile:
        outfile.write(citation)

def PrintCitation(d=None):
    if d is not None: WriteCitation(d)
    print()
    print(print_citation)  

def PrintUnderline(text, qHeavy=False):
    print(("\n" + text))
    n = len(text)
    if text.startswith("\n"): n -= 1
    print((("=" if qHeavy else "-") * n))

def FlowText(text, n=60):
    """Split text onto lines of no more that n characters long
    """
    lines = ""
    while len(text) > 0:
        if len(lines) > 0: lines += "\n"
        if len(text) > n:
            # split at no more than 60
            iEnd = n
            while iEnd > 0 and text[iEnd] != " ": iEnd-=1
            if iEnd == 0:
                # there was nowhere to split it at a blank, just have to split at 60
                lines += text[:n]
                text = text[n:]
            else:
                # split at blank
                lines += text[:iEnd]
                text = text[iEnd+1:]  # skip blank
        else:
            lines += text
            text = ""
    return lines
        
class Finalise(object):
    def __enter__(self):
        pass
    def __exit__(self, type, value, traceback):
        ptm = parallel_task_manager.ParallelTaskManager_singleton()
        ptm.Stop()

def writerow(fh, row):
    # CSV format specifies CRLF line endings: https://tools.ietf.org/html/rfc4180
    fh.write("\t".join(map(str, row)) + "\r\n")

def getrow(row):
    # CSV format specifies CRLF line endings: https://tools.ietf.org/html/rfc4180
    return "\t".join(map(str, row)) + "\r\n"

def version_parse_simple(sem_version):
    return list(map(int, sem_version.split(".")[:3]))

def file_open(filename_with_gz, mode, gz):
    if gz:
        return gzip.open(filename_with_gz + ".gz", mode)
    else:
        return open(filename_with_gz, mode)
