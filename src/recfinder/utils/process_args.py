import os
from typing import Optional
from . import util
import shutil
from .. import nThreadsDefault, __version__
from .. import helpinfo

class Options(object):#
    def __init__(self):
        
        self.tree_program = "iqtree"
        self.sequence_type = "CODON11"
        self.evolution_model = "TEST"
        self.qStartFromMSA = False
        self.nProcessAlg = None
        self.nMSATree = nThreadsDefault
        self.name = None   # name to identify this set of results
        self.extended_filename = False
        self.partial_constraints = None
        self.output_prefix = None
        self.gene_tree = None
        self.root_node = None
        self.nalign = "1000"
        self.outgroups = None 
        self.species_of_interest = None
        self.usr_node_aln = None 
        self.usr_state = None 
        self.usr_tree = None

        

    def what(self):
        for k, v in self.__dict__.items():
            if v == True:
                print(k)

def GetDirectoryArgument(arg: str) -> str:
    directory = os.path.abspath(arg)
    if not os.path.isfile(directory) and directory[-1] != os.sep: 
        directory += os.sep
    if not os.path.exists(directory):
        print("Specified directory doesn't exist: %s" % directory)
        util.Fail()
    return directory


def GetFileArgument(arg: str) -> str:
    file_path = os.path.abspath(arg)
    if not os.path.isfile(file_path):
        print("Specified file doesn't exist: %s" % file_path)
        util.Fail()
    return file_path

def ProcessArgs(args: list[str]):

    if len(args) == 0 or args[0] == "--help" or args[0] == "help" or args[0] == "-h":
        helpinfo.PrintHelp()
        util.Success()

    if args[0] == "-v" or args[0] == "--version":
        print(f"RecFinder:v{__version__}")
        util.Success()

    options = Options()
    alnDir = None
    resultsDir_nonDefault = None
    q_selected_tree_options = False
    
    """
    -f: store alnDir
    -b: store workingDir
    -fg: store orthologuesDir 
    -ft: store orthologuesDir 
    + xml: speciesXMLInfoFN
    """    
    
    while len(args) > 0:
        arg = args.pop(0)

        if arg == "-f" or arg == "--faln":
            if options.qStartFromMSA:
                print("Repeated argument: -f/--faln\n")
                util.Fail()

            options.qStartFromMSA = True

            if len(args) == 0:
                print("Missing option for command line argument %s" % arg)
                util.Fail()
            alnDir = GetDirectoryArgument(args.pop(0))

        elif arg == "-una" or arg == "--usr-node-aln":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()

            options.usr_node_aln = GetFileArgument(args.pop(0))

        
        elif arg == "-us" or arg == "--usr-state":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()

            options.usr_state = GetFileArgument(args.pop(0))

        
        elif arg == "-ut" or arg == "--usr-tree":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()

            options.usr_tree = GetFileArgument(args.pop(0))


        elif arg == "--outgroups":
            arg_outgroup = arg
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            isfile = False
            arg = args.pop(0)
            try:
                if os.path.isfile(arg):
                    isfile = True
                    outgroup_path = GetFileArgument(arg)
                    with open(outgroup_path) as reader:
                        outgroups = []
                        for line in reader:
                            line = line.replace("\n", "").strip()
                            outgroups.append(line)
                else:
                    if "," in arg:
                        outgroups = [og.strip() for og in arg.split(",")]
                    else:
                        outgroups = [arg]
            except:
                print("Invalid argument for option %s: %s" % (arg_outgroup, arg))
                print('Valid options are for instance "Nuphar japonica,Nymphaea jamesoniana"\n')
                util.Fail()

            if len(outgroups) == 0 and not isfile:
                print("Missing option for command line argument from file %s\n" % arg)
                util.Fail()
                        
            options.outgroups = outgroups

        elif arg == "--spofin":
            arg_species = arg
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            isfile = False
            arg = args.pop(0)
            try:
                if os.path.isfile(arg):
                    spofin_path = GetFileArgument(arg)
                    isfile = True
                    with open(spofin_path) as reader:
                        species_of_interest = []
                        for line in reader:
                            line = line.replace("\n", "").strip()
                            species_of_interest.append(line)
                else:
                    if "," in arg:
                        species_of_interest = [og.strip() for og in arg.split(",")]
                    else:
                        species_of_interest = [arg]
            except:
                print("Invalid argument for option %s: %s" % (species, arg))
                print('Valid options are for instance "Talipariti hamabo,Aegilops geniculata"\n')
                util.Fail()
                
            if len(species_of_interest) == 0 and not isfile:
                print("Missing option for command line argument from file %s\n" % arg)
                util.Fail()
                        
            options.species_of_interest = species_of_interest

        elif arg == "-st":
            arg_st = arg
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()

            arg = args.pop(0).upper()
            try:
                if "CODON" in arg:
                    arg1, arg2 = arg[:5], int(arg[5:])
                    if arg1 == "CODON" and arg2 in [*range(1, 26)]:
                        options.sequence_type = arg
                elif "AA" in arg:
                    options.sequence_type = arg
            except:
                print("Invalid argument for option %s: %s" % (arg_st, arg))
                print("Valid options are for instance 'CODON[1-11]' or 'AA'\n")
                print("For more information please refer to http://www.iqtree.org/doc/Substitution-Models#codon-models")
                util.Fail()

        elif arg == "-m" or arg == "--evolution-model":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            options.evolution_model = args.pop(0).upper()

        elif arg == "-te":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            options.gene_tree = args.pop(0)

        elif arg == "-g":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            options.partial_constraints = args.pop(0)

        elif arg == "--alisim":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()   
            options.output_prefix = args.pop(0)
        
        elif arg == "-m":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()   
            options.evolution_model = args.pop(0)  

        elif arg == "--root-seq":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()   
            options.root_node = args.pop(0)

        elif arg == "--num-alignments":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()   
            options.root_node = str(int(args.pop(0)))          

        elif arg == "-t" or arg == "--threads":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)

            try:
                options.nMSATree = int(arg)
            except:
                print("Incorrect argument for number of IQTREE threads: %s\n" % arg)
                util.Fail()

        elif arg == "-a" or arg == "--algthreads":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            try:
                options.nProcessAlg = int(arg)
            except:
                print("Incorrect argument for number of RecFinder Analysis threads: %s\n" % arg)
                util.Fail()  

        elif arg == "--save-space":
            options.save_space = True

        elif arg == "-n" or arg == "--name":  
            if options.name:
                print("Repeated argument: -n/--name")
                util.Fail()
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail() 
            options.name = args.pop(0)
            while options.name.endswith("/"): options.name = options.name[:-1]
            if any([symbol in options.name for symbol in [" ", "/"]]): 
                print("Invalid symbol for command line argument %s\n" % arg)
                util.Fail()

        elif arg == "-o" or arg == "--output":  
            if resultsDir_nonDefault != None:
                print("Repeated argument: -o/--output")
                util.Fail()
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            resultsDir_nonDefault = args.pop(0)

            while resultsDir_nonDefault.endswith("/"): 
                resultsDir_nonDefault = resultsDir_nonDefault[:-1]

            resultsDir_nonDefault += "/"
            if os.path.exists(resultsDir_nonDefault):
                print("ERROR: non-default output directory already exists: %s\n" % resultsDir_nonDefault)
                util.Fail()

            if " " in resultsDir_nonDefault:
                print("ERROR: non-default output directory cannot include spaces: %s\n" % resultsDir_nonDefault)
                util.Fail()
            checkDirName = resultsDir_nonDefault

            while checkDirName.endswith("/"):
                checkDirName = checkDirName[:-1]

            path, newDir = os.path.split(checkDirName)
            if path != "" and not os.path.exists(path):
                print("ERROR: location '%s' for results directory '%s' does not exist.\n" % (path, newDir))
                util.Fail()

        elif arg == "-h" or arg == "--help":
            helpinfo.PrintHelp()
            util.Success()

        elif arg == "-efn" or arg == "--extended-filename":
            options.extended_filename = True

        else:
            print("Unrecognised argument: %s\n" % arg)
            util.Fail()

    # set a default for number of algorithm threads
    if options.nProcessAlg is None:
        options.nProcessAlg = min(16, max(1, int(options.nMSATree/8)))

    if options.nMSATree < 1:
        print("ERROR: Number of '-t' threads cannot be fewer than 1, got %d" % options.nMSATree)
        util.Fail()    

    if options.nProcessAlg < 1:
        print("ERROR: Number of '-a' threads cannot be fewer than 1, got %d" % options.nProcessAlg)
        util.Fail()  

    # check argument combinations       
    if not options.qStartFromMSA :
        print("ERROR: Please specify the input directory for RecFinder using the option: '-f'.")
        util.Fail()

    print()
    util.PrintTime("Starting RecFinder v%s" % __version__)    
    print("%d thread(s) for highly parallel tasks (IQ-TREE inference etc.)" % options.nMSATree)
    print("%d thread(s) for RecFinder algorithm\n" % options.nProcessAlg)

    return options, alnDir, resultsDir_nonDefault


def DeleteDirectoryTree(d):
    if os.path.exists(d): 
        try:
            shutil.rmtree(d)
        except OSError:
            time.sleep(1)
            shutil.rmtree(d, True)   


def CheckOptions(options, speciesToUse):
    """Check any optional arguments are valid once we know what species are in the analysis
    - user supplied species tree
    """
    if options.speciesTreeFN:
        expSpecies = list(SpeciesNameDict(files.FileHandler.GetSpeciesIDsFN()).values())
        orthologues.CheckUserSpeciesTree(options.speciesTreeFN, expSpecies)

    # check can open enough files
    n_extra = 50
    q_do_orthologs = not any((options.qStopAfterPrepare, options.qStopAfterGroups, options.qStopAfterSeqs, options.qStopAfterAlignments, options.qStopAfterTrees))
    if q_do_orthologs:
        n_sp = len(speciesToUse)
        wd = files.FileHandler.GetWorkingDirectory_Write()
        wd_files_test = wd + "Files_test/"
        fh = []
        try:
            if not os.path.exists(wd_files_test):
                os.mkdir(wd_files_test)
            for i_sp in range(n_sp):
                di = wd_files_test + "Sp%d/" % i_sp
                if not os.path.exists(di):
                    os.mkdir(di)
                for j_sp in range(1):  # We only create a linear number of ortholog files now
                    fnij = di + "Sp%d.txt" % j_sp
                    fh.append(open(fnij, 'w'))
            # create a few extra files to be safe
            for i_extra in range(n_extra):
                fh.append(open(wd_files_test + "Extra%d.txt" % i_extra, 'w'))
            # close the files again and delete
            for fhh in fh:
                fhh.close()
            DeleteDirectoryTree(wd_files_test)
        except IOError as e:
            if str(e).startswith("[Errno 24] Too many open files"):
                util.number_open_files_exception_advice(len(speciesToUse), False)
                for fhh in fh:
                    fhh.close()
                DeleteDirectoryTree(wd_files_test)
                util.Fail()
            else:
                for fhh in fh:
                    fhh.close()
                DeleteDirectoryTree(wd_files_test)
                print("ERROR: Attempted to open required files for RecFinder run but an unexpected error occurred. \n\nStacktrace:")
                raise
    return options