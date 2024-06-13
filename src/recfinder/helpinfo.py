
from . import nThreadsDefault
from .citation import print_citation

def PrintHelp():  

    print("")
    print("SIMPLE USAGE:") 
    print("Run full OrthoFinder analysis on FASTA format proteomes in <dir>")
    print("  orthofinder [options] -f <dir>")   
    # print("")
    # print("Add new species in <dir1> to previous run in <dir2> and run new analysis")
    # print("  orthofinder [options] -f <dir1> -b <dir2>")
    print("")
    print("To assign species from <dir1> to existing OrthoFinder orthogroups in <dir2>")
    print("  orthofinder [options] --assign <dir1> --core <dir2>")
    print("") 
      
    print("OPTIONS:")
    print(" -t <int>                Number of parallel sequence search threads [Default = %d]" % nThreadsDefault)
    print(" -a <int>                Number of parallel analysis threads")
    print(" -te <file>              Fixing user tree. That means, no tree search is performed and IQ-TREE computes the log-likelihood of the fixed user tree.")
    print(" --save-space            Only create one compressed orthologs file per species")
    print(" -g <file>               Specify a topological constraint tree file in NEWICK format. The constraint tree can be a multifurcating tree and need not to include all taxa.")
    print(" -p <dir>                Write the temporary pickle files to <dir>")
    print(" --num-alignments <int>       Set the number of output datasets. Default: 1000")
    print(" -m <MODEL>              Specify the model name. See Substitution Models (http://www.iqtree.org/doc/Substitution-Models) for the list of supported models.")
    print(" --alisim <file??>	       Activate AliSim and specify the output alignment filename.")
    print(" -z                      Don't trim MSAs (columns>=90% gap, min. alignment length 500)")
    print(" -n <txt>                Name to append to the results directory")  
    print(" -o <txt>                Non-default results directory")  
    print(" -h                      Print this help text")
    print(" -efn                    Extend the output directory name with the name of the codon model")


    print("")    
    print("WORKFLOW STOPPING OPTIONS:")   
    print(" -op                     Stop after preparing input files for BLAST" )
    print(" -og                     Stop after inferring orthogroups")
    print(" -os                     Stop after writing sequence files for orthogroups")
    print("                         (requires '-M msa')")
    print(" -oa                     Stop after inferring alignments for orthogroups")
    print("                         (requires '-M msa')")
    print(" -ot                     Stop after inferring gene trees for orthogroups " )
   
    print("")   
    print("WORKFLOW RESTART COMMANDS:") 
    print(" -b  <dir>               Start OrthoFinder from pre-computed BLAST results in <dir>")   
    print(" -fg <dir>               Start OrthoFinder from pre-computed orthogroups in <dir>")
    print(" -ft <dir>               Start OrthoFinder from pre-computed gene trees in <dir>")

    print("")
    print("VERSION:")
    print(" -v                      Show the current version number")
    
    print("")
    print("LICENSE:")
    print(" Distributed under the GNU General Public License (GPLv3). See License.md")
    print(print_citation) 
