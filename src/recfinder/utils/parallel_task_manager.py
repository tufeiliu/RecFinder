# -*- coding: utf-8 -*-

import os
import sys
import time
import types
import datetime
import traceback
import subprocess
import contextlib
import threading
import concurrent.futures
import multiprocessing as mp
from . import util, files
from concurrent.futures import ProcessPoolExecutor, wait, as_completed
from .. import my_env
try:
    import queue
except ImportError:
    import Queue as queue     
from typing import Optional, Union, Callable, Any, List, Dict, Tuple
# uncomment to get round problem with python multiprocessing library that can set all cpu affinities to a single cpu
# This can cause use of only a limited number of cpus in other cases so it has been commented out
# if sys.platform.startswith("linux"):
#     with open(os.devnull, "w") as f:
#         subprocess.call("taskset -p 0xffffffffffff %d" % os.getpid(), shell=True, stdout=f)

q_print_first_traceback_0 = False
lock = threading.Lock()


PY2 = sys.version_info <= (3,)

def print_traceback(e):
    if PY2:
        traceback.print_exc()
    else:
        traceback.print_tb(e.__traceback__)


def PrintTime(message):
    print((str(datetime.datetime.now()).rsplit(".", 1)[0] + " : " + message))
    sys.stdout.flush()

def PrintNoNewLine(text):
    sys.stdout.write(text)

def ParallelComputer(func: Callable[..., Any], task_queue: mp.Queue, nthreads: int, *args):
    
    result_queue = mp.Queue()
    queue_len = task_queue.qsize()
    
    for _ in range(nthreads):
        task_queue.put(None)

    processes = [
        mp.Process(target=func, args=(task_queue, result_queue, *args))
        for i_ in range(nthreads)]
    for proc in processes:
        proc.start()

    results = [result_queue.get() for _ in range(queue_len)]
    ManageQueue(processes, task_queue)
    
    return results


def ParallelWriter(func, task_queue, nthreads, *args):
    
    result_queue = mp.Queue()
    mplock = mp.Lock()

    for _ in range(nthreads):
        task_queue.put(None)
    
    processes = []      
    for _ in range(nthreads):
        worker_process = mp.Process(target=func, args=(task_queue, result_queue, mplock, *args))
        processes.append(worker_process)
        worker_process.start()

    while not result_queue.empty():
        print(result_queue.get())

    for worker_process in processes:    
        worker_process.join()

    ManageQueue(processes, task_queue)

def ParellelReader(func, task_queue, nthreads, *args):

    result_queue = mp.Queue()
    queue_len = task_queue.qsize()
    
    for _ in range(nthreads):
        task_queue.put(None)

    with concurrent.futures.ThreadPoolExecutor(max_workers=nthreads) as executor:
        futures = [executor.submit(func, task_queue, result_queue) for _ in range(nthreads)]
    
    concurrent.futures.wait(futures)

    return result_queue


def ManageQueue(runningProcesses, cmd_queue):
    """Manage a set of runningProcesses working through cmd_queue.
    If there is an error the exit all processes as quickly as possible and 
    exit via Fail() methods. Otherwise return when all work is complete
    """            
    # set all completed processes to None
    qError = False
#    dones = [False for _ in runningProcesses]
    nProcesses = len(runningProcesses)
    nProcesses_list = list(range(nProcesses))
    while True:
        if runningProcesses.count(None) == len(runningProcesses): break
        time.sleep(.1)
#        for proc in runningProcesses:
        for i in nProcesses_list:
            proc = runningProcesses[i]
            if proc == None: continue
            if not proc.is_alive():
                if proc.exitcode != 0:
                    qError = True
                    while True:
                        try:
                            cmd_queue.get(True, .1)
                        except queue.Empty:
                            break
                runningProcesses[i] = None
    if qError:
        Fail()


def RunCommand_Simple(command):
    subprocess.call(command, env=my_env, shell=True)


def CanRunCommand(command, qAllowStderr = False, qPrint = True, qRequireStdout=True, qCheckReturnCode=False):
    if qPrint: PrintNoNewLine("Test can run \"%s\"" % command)       # print without newline
    capture = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=my_env)
    capture.wait()
    stdout = [x for x in capture.stdout]
    stderr = [x for x in capture.stderr]
    if qCheckReturnCode:
        return_code_check = (capture.returncode == 0)
    else:
        return_code_check = True
    if (len(stdout) > 0 or not qRequireStdout) and (qAllowStderr or len(stderr) == 0) and return_code_check:
        if qPrint: print(" - ok")
        return True
    else:
        if qPrint: print(" - failed")
        if not return_code_check:
            print("Returned a non-zero code: %d" % capture.returncode)
        print("\nstdout:")        
        for l in stdout: print(l)
        print("\nstderr:")        
        for l in stderr: print(l)
        return False



q_print_first_traceback_1 = False
def Worker_RunMethod(Function, args_queue):
    while True:
        try:
            args = args_queue.get(True, .1)
            Function(*args)
        except queue.Empty:
            return 
        except Exception as e:
            with lock:
                print("Error in function: " + str(Function))
                global q_print_first_traceback_1
                if not q_print_first_traceback_1:
                    print_traceback(e)
                    q_print_first_traceback_1 = True
            return
        except:
            with lock:
                print("WARNING: Unknown caught unknown exception")

def RunMethodParallel(Function, args_queue, nProcesses):
    runningProcesses = [mp.Process(target=Worker_RunMethod, args=(Function, args_queue)) for i_ in range(nProcesses)]
    for proc in runningProcesses:
        proc.start()
    ManageQueue(runningProcesses, args_queue)


def _I_Spawn_Processes(message_to_spawner, message_to_PTM):
    """
    Args:
        message_queue - for passing messages that a new queue of tasks should be started (PTM -> I_Space_Processes) or that the tasks are complete
        cmds_queue - queue containing tasks that should be done
    Use:
        A process should be started as early as possible (while RAM usage is low) with this method as its target.
        This is now a separate process with low RAM usage.
        Each time some parallel work is required then the queue for that is placed in the message_queue by the PTM.
        _I_Spawn_Processes - will spawn parallel processes when instructed by the message_queue in the message_queue and get them
        working on the queue. When the queue is empty it will wait for the next one. It can receive a special signal to exit - the None
        object
    """
    while True:
        try:
            # peak in qoq - it is the only method that tried to remove things from the queue
            message = message_to_spawner.get(timeout=.1)
            if message is None:
                # Respond to request to terminate
                return
            # In which case, thread has been informed that there are tasks in the queue.
            func, args_list, n_parallel = message
            futures = []
            n_to_do = len(args_list)
            # Version 1: n worker threads for executing the method and a list of N arguments for calling the method
            with ProcessPoolExecutor(n_parallel) as pool:
                for args in args_list:
                    futures.append(pool.submit(func, *args))
                # for i, _ in as_completed(futures):
                #     n_done = i+1
                #     if n_done >= 0 and divmod(n_done, 10 if n_done <= 200 else 100 if n_done <= 2000 else 1000)[1] == 0:
                #         PrintTime("Done %d of %d" % (n_done, n_to_do))
            # Version 2: launch n worker threads each executing a worker method that takes tasks from a queue
            with ProcessPoolExecutor(n_parallel) as pool:
                for args in args_list:
                    futures.append(pool.submit(func, *args))
            wait(futures)
            message_to_PTM.put("Done")
            time.sleep(1)

        except KeyboardInterrupt:
            # Handle KeyboardInterrupt to exit gracefully
            message_to_PTM.put("KeyboardInterrupt")
            return

        except queue.Empty:
            time.sleep(1)  # there wasn't anything this time, sleep then try again
    pass

class ParallelTaskManager_singleton:
    """
    Creating new process requires forking parent process and can lea to very high RAM usage. One way to mitigate this is
    to create the pool of processes as early in execution as possible so that the memory footprint is low. The
    ParallelTaskManager takes care of that, and can be used by calling `RunParallelOrderedCommandLists` above.
    Apr 2023 Update:
    When running external programs there is no need to use multiprocessing, multithreading is sufficient since new process
    will be created anyway, so the SIL is no longer an issue.
    """
    class __Singleton(object):
        def __init__(self):
            """Implementation:
            Allocate a thread that will perform all the tasks
            Communicate with it using a queue.
            When provided with a list of commands it should fire up some workers and get them to run the commands and then exit.
            An alternative would be they should always stay alive - but then they could die for some reason? And I'd have to check how many there are.
            """
            self.message_to_spawner = mp.Queue()
            self.message_to_PTM = mp.Queue()
            # Orders/Messages:
            # None (PTM -> spawn_thread) - thread should return (i.e. exit)
            # 'Done' (spawn_thread -> PTM) - the cmds from the cmd queue have completed
            # Anything else = (nParallel, nTasks) (PTM -> spawn_thread) - cmds (nTasks of them) have been placed in the cmd queue,
            #   they should be executed using nParallel threads
            self.manager_process = mp.Process(target=_I_Spawn_Processes, args=(self.message_to_spawner, self.message_to_PTM))
            self.manager_process.start()
    instance = None

    def __init__(self):
        if not ParallelTaskManager_singleton.instance:
            ParallelTaskManager_singleton.instance = ParallelTaskManager_singleton.__Singleton()

    def RunParallel(self, func, args_list, nParallel):
        """
        Args:
            cmd_list - list of commands or list of lists of commands (in which elements in inner list must be run in order)
            nParallel - number of parallel threads to use
            qShell - should the tasks be run in a shell
        """
        self.instance.message_to_spawner.put((func, args_list, nParallel))
        while True:
            try:
                signal = self.instance.message_to_PTM.get()
                if signal == "Done":
                    return
            except queue.Empty:
                pass
            time.sleep(1)

    def Stop(self):
        """Warning, cannot be restarted"""
        self.instance.message_to_spawner.put(None)
        self.instance.manager_process.join()

# ================================================================================
def RunParallelMethods(func, args_list, nProcesses):
    """nProcesss - the number of processes to run in parallel
    commands - list of lists of commands where the commands in the inner list are completed in order (the i_th won't run until
    the i-1_th has finished).
    """
    ptm = ParallelTaskManager_singleton()
    ptm.RunParallel(func, args_list, nProcesses)

def Success():
    ptm = ParallelTaskManager_singleton()
    ptm.Stop()
    sys.exit()

def Fail():
    sys.stderr.flush()
    ptm = ParallelTaskManager_singleton()
    ptm.Stop()
    print("ERROR: An error occurred, ***please review the error messages*** they may contain useful information about the problem.")
    sys.exit(1)



def RunCommand(command, qPrintOnError=False, qPrintStderr=True):

    """ Run a single command """

    # popen = subprocess.Popen(command, env=my_env, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # if qPrintOnError:
    #     stdout, stderr = popen.communicate()
    #     if popen.returncode != 0:
    #         print(("\nERROR: external program called by RecFinder returned an error code: %d" % popen.returncode))
    #         print(("\nCommand: %s" % command))
    #         print(("\nstdout:\n%s" % stdout))
    #         print(("stderr:\n%s" % stderr))
    #     elif qPrintStderr and len(stderr) > 0 and not stderr_exempt(stderr):
    #         print("\nWARNING: program called by RecFinder produced output to stderr")
    #         print(("\nCommand: %s" % command))
    #         print(("\nstdout:\n%s" % stdout))
    #         print(("stderr:\n%s" % stderr))
    #     return popen.returncode
    # else:
    #     popen.communicate()
    #     return popen.returncode


    try:
        popen = subprocess.Popen(command, env=my_env, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = popen.communicate()

        if qPrintOnError and popen.returncode != 0:
            with lock:
                print(f"\nERROR: external program called by RecFinder returned an error code: {popen.returncode}")
                print(f"\nCommand: {command}")
                print(f"\nstdout:\n{stdout.decode()}")
                print(f"stderr:\n{stderr.decode()}")

        elif qPrintStderr and len(stderr) > 0 and not stderr_exempt(stderr.decode()):
            with lock:
                print("\nWARNING: program called by RecFinder produced output to stderr")
                print(f"\nCommand: {command}")
                print(f"\nstdout:\n{stdout.decode()}")
                print(f"stderr:\n{stderr.decode()}")

        return popen.returncode

    except Exception as e:
        print(f"Exception occurred while running command: {command}")
        print(f"Exception: {e}")
        return -1

def Worker_RunCommands_And_Move(cmd_and_filename_queue, 
                                nProcesses, nToDo, 
                                qListOfLists, q_print_on_error, 
                                q_always_print_stderr,
                                fileDir, 
                                convert_statefiles,
                                sequence_type, 
                                delete_files,
                                files_to_keep,
                                files_to_remove):
    """
    Continuously takes commands that need to be run from the cmd_and_filename_queue until the queue is empty. If required, moves 
    the output filename produced by the cmd to a specified filename. The elements of the queue can be single cmd_filename tuples
    or an ordered list of tuples that must be run in the provided order.
  
    Args:
        cmd_and_filename_queue - queue containing (cmd, actual_target_fn) tuples (if qListOfLists is False) or a list of such 
            tuples (if qListOfLists is True). Alternatively, 'cmd' can be a python fn and actual_target_fn the fn to call it on.
        nProcesses - the number of processes that are working on the queue.
        nToDo - The total number of elements in the original queue
        qListOfLists - Boolean, whether each element of the queue corresponds to a single command or a list of ordered commands
        qShell - Boolean, should a shell be used to run the command.
        
    Implementation:
        nProcesses and nToDo are used to print out the progress.
    """
    while True:
        try:
            i, command_fns_list = cmd_and_filename_queue.get(True, 1)
            nDone = i - nProcesses + 1
            if nDone >= 0 and divmod(nDone, 10 if nToDo <= 200 else 100 if nToDo <= 2000 else 1000)[1] == 0:
                PrintTime("Done %d of %d" % (nDone, nToDo))
            if not qListOfLists:
                command_fns_list = [command_fns_list]
            for command, fns in command_fns_list:
                if isinstance(command, types.FunctionType):
                    # This will block the process, but it is ok for trimming, it takes minimal time
                    fn = command
                    fn(fns)
                else:
                    if not isinstance(command, str):
                        with lock:
                            print("ERROR: Cannot run command: " + str(command))
                            print("Please report this issue.")
                    else:
                        RunCommand(command, qPrintOnError=q_print_on_error, qPrintStderr=q_always_print_stderr)
                        if fns != None:
                            actual, target = fns
                            if os.path.exists(actual):
                                wait_for_file_completion(actual)
                                os.rename(actual, target)

            cmd_and_filename_queue.task_done()

            if convert_statefiles:
                handle_state_files(fileDir, sequence_type, delete_files, command)
            if delete_files:
                clean_up_files(fileDir, files_to_keep, files_to_remove)
                
        except queue.Empty:
            return  
        except KeyboardInterrupt:
            with lock:
                print("Process interrupted by user. Cleaning up...")
            break   
        except Exception as e:
            with lock:
                print("WARNING: ")
                print(str(e))
                global q_print_first_traceback_0
                if not q_print_first_traceback_0:
                    print_traceback(e)
                    q_print_first_traceback_0 = True
        except:
            with lock:
                print("WARNING: Unknown caught unknown exception")


def RunParallelCommands(nProcesses, commands, fileDir, convert_statefiles=False, sequence_type=None,
                        delete_files=False,  files_to_keep = None, files_to_remove = None,
                        qListOfList=False, q_print_on_error=False, q_always_print_stderr=False):
    if qListOfList:
        commands_and_no_filenames = [[(cmd, None) for cmd in cmd_list] for cmd_list in commands]
    else:
        commands_and_no_filenames = [(cmd, None) for cmd in commands]
    RunParallelCommandsAndMoveResultsFile(nProcesses, commands_and_no_filenames,
                                          fileDir, convert_statefiles, sequence_type, delete_files, files_to_keep, files_to_remove,
                                          qListOfList, q_print_on_error,
                                          q_always_print_stderr)


def RunParallelCommandsAndMoveResultsFile(nProcesses, 
                                          commands_and_filenames, 
                                          fileDir,
                                          convert_statefiles=False,
                                          sequence_type = None, 
                                          delete_files=False,
                                          files_to_keep=None,
                                          files_to_remove=None,
                                          qListOfList=False, 
                                          q_print_on_error=False,
                                          q_always_print_stderr=False):
    """
    Calls the commands in parallel and if required moves the results file to the required new filename
    Args:
        nProcess - the number of parallel process to use
        commands_and_filenames : tuple (cmd, actual_target_fn) where actual_target_fn = None if no move is required 
                                 and actual_target_fn = (actual_fn, target_fn) is actual_fn is produced by cmd and this 
                                 file should be moved to target_fn
        actual_target_fn - None if the cmd will save the results file to outfilename_proposed 
                           otherwise (actual_fn, outfilename_proposed)
        qListOfList - if False then commands_and_filenames is a list of (cmd, actual_target_fn) tuples
                      if True then commands_and_filenames is a list of lists of (cmd, actual_target_fn) tuples where the elements 
                      of the inner list need to be run in the order they appear.
        q_print_on_error - If error code returend print stdout & stederr
    """
    cmd_queue = queue.Queue()
    i = -1
    for i, cmd in enumerate(commands_and_filenames):
        cmd_queue.put((i, cmd))
    try:
        with concurrent.futures.ThreadPoolExecutor() as executor:
            futures = [executor.submit(Worker_RunCommands_And_Move,
                                    cmd_queue,
                                    nProcesses,
                                    i+1,
                                    qListOfList,
                                    q_print_on_error,
                                    q_always_print_stderr,
                                    fileDir,
                                    convert_statefiles=convert_statefiles,
                                    sequence_type=sequence_type, 
                                    delete_files=delete_files,
                                    files_to_keep=files_to_keep,
                                    files_to_remove=files_to_remove)
                    for _ in range(nProcesses)]

        cmd_queue.join()
        concurrent.futures.wait(futures)

    except KeyboardInterrupt:
        with lock:
            print("Process interrupted by user. Cleaning up...")
        for future in futures:
            future.cancel()
        cmd_queue.queue.clear()



def wait_for_file_completion(filepath, wait_time=1, retries=10):
    """Wait for a file to be fully written by checking if its size remains constant."""
    for attempt in range(retries):
        if os.path.exists(filepath):
            previous_size = -1
            for sub_attempt in range(retries):
                current_size = os.path.getsize(filepath)
                if current_size == previous_size:
                    # print(f"File {filepath} is ready after {attempt * wait_time + sub_attempt * wait_time} seconds")
                    return True
                previous_size = current_size
                time.sleep(wait_time)
        else:
            print(f"File {filepath} does not exist, retrying... ({attempt + 1}/{retries})")
            time.sleep(wait_time)
    raise FileNotFoundError(f"File {filepath} not found after {wait_time * retries} seconds.")

def handle_state_files(fileDir, sequence_type, delete_files, command):
    combined_prot_seqs_d = os.path.join(os.path.dirname(os.path.dirname(fileDir)), "Alignment_Files")
    
    if not os.path.exists(combined_prot_seqs_d):
        os.mkdir(combined_prot_seqs_d)
    
    phy_path = command.split()[2]
    # print(f"Phy path: {phy_path}")

    if not os.path.isfile(phy_path):
        print(f"ERROR: Phy file not found: {phy_path}")
        return

    processed_files = set()

    for file in os.listdir(fileDir):
        file_path = os.path.join(fileDir, file)
        # print(f"Processing file: {file_path}")

        if os.path.isfile(file_path) and file_path not in processed_files:
            file_extension = file.rsplit(".", 1)[1].lower()
            if file_extension == "state":
                try:
                    wait_for_file_completion(file_path)
                    
                    # print(f"Reading phy file: {phy_path}")
                    alignment_dict = files.FileReader.ReadPhyFile(phy_path)
                    
                    # print(f"Reading state file: {file_path}")
                    node_seq_dict = files.FileReader.ReadStateFile(file_path)
                    
                    combined_seq_dict = {k: v for d in (node_seq_dict, alignment_dict) for k, v in d.items()}
                    
                    if sequence_type == "AA":
                        combined_prot_seqs_dict = combined_seq_dict
                    else:
                        combined_prot_seqs_dict, _ = util.GetSeqsDict(combined_seq_dict, sequence_type)
                    
                    combined_prot_seqs_fn = os.path.join(combined_prot_seqs_d, file.rsplit(".")[0] + ".aln")
                    # print(f"Writing combined sequences to: {combined_prot_seqs_fn}")
                    files.FileWriter.WriteSeqsToAln(combined_prot_seqs_dict, combined_prot_seqs_fn)
                    
                    processed_files.add(file_path)

                    if delete_files:
                        with contextlib.suppress(FileNotFoundError):
                            os.remove(file_path)
                            print(f"Deleted file: {file_path}")

                except FileNotFoundError as e:
                    print(f"ERROR: File not found: {e.filename}")
                except Exception as e:
                    print(f"ERROR processing file {file_path}: {e}")
    # return processed_files


def clean_up_files(fileDir, files_to_keep, files_to_remove):
    if files_to_keep:
        for file in os.listdir(fileDir):
            file_path = os.path.join(fileDir, file)
            if os.path.isfile(file_path):
                if file.rsplit(".", 1)[-1].lower() not in files_to_keep:
                    with contextlib.suppress(FileNotFoundError):
                        os.remove(file_path)

    elif files_to_remove:
        for file in os.listdir(fileDir):
            file_path = os.path.join(fileDir, file)
            if os.path.isfile(file_path):
                if file.rsplit(".", 1)[-1].lower() in files_to_remove:
                    with contextlib.suppress(FileNotFoundError):
                        os.remove(file_path)




