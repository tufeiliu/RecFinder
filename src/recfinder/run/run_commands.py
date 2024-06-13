from ..utils import util, files, parallel_task_manager
import shutil
from .. import my_env
import os
from typing import Optional, Dict, List, Tuple


build_gene_tree_command_str = "iqtree2 -s alignment_file -T nthreads -asr -m evolution_model -pre path_to_output" 
Monte_Carlo_simulation_command_str =  "iqtree2 --alisim output_prefix -T nthreads -m best_evolution_model -te gene_tree --root-seq root_node --num-alignments nalign"


def GetGeneTreeBuildCommand(alignment_file: str, 
                            nthreads: int, 
                            evolution_model: str, 
                            path_to_output: str,
                            sequence_type: Optional[str] = None,
                            gene_tree: Optional[str] = None, 
                            partial_constraints: Optional[str] = None) -> str:

    command = build_gene_tree_command_str.\
                                    replace("alignment_file", alignment_file).\
                                    replace("nthreads", str(nthreads)).\
                                    replace("evolution_model", evolution_model).\
                                    replace("path_to_output", path_to_output).split()

    if sequence_type:
        command.extend(["-st", sequence_type])
    if gene_tree:
        command.extend(["-te", gene_tree])
    if partial_constraints:
        command.extend(["-g", partial_constraints])
    
    return " ".join(command)


def GetMCsimulationCommand(output_prefix: str, 
                           nthreads: int,
                           best_evolution_model: str, 
                           gene_tree: str, 
                           root_node: str, 
                           nalign: str) -> str:
    command = Monte_Carlo_simulation_command_str.\
                                     replace("output_prefix", output_prefix).\
                                     replace("nthreads", str(nthreads)).\
                                     replace("best_evolution_model", best_evolution_model).\
                                     replace("gene_tree", gene_tree).\
                                     replace("root_node", root_node).\
                                     replace("nalign", nalign)
    
    return [command]


 
def RunCommand(nthreads: int, commands: list[str], 
               print_info: str, fileDir: str, 
               convert_statefiles: bool = False,
               sequence_type: Optional[str] = None,
               delete_files: bool = False,
               files_to_keep: Optional[List[str]] = None, 
               files_to_remove: Optional[List[str]] = None) -> None:

    parallel_task_manager.RunParallelCommands(nthreads, commands, fileDir,
                                              convert_statefiles = convert_statefiles,
                                              sequence_type = sequence_type,
                                              delete_files = delete_files,
                                              files_to_keep = files_to_keep,
                                              files_to_remove = files_to_remove, 
                                              qListOfList = False,
                                              q_print_on_error = False, 
                                              q_always_print_stderr = False)
    util.PrintTime(print_info)
    print("Using %d thread(s)" % nthreads)
    util.PrintTime("This may take some time....")