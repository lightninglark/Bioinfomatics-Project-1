'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
from setuptools.dist import sequence
Created on Apr 21, 2017

Project 1: Project 1 ADT is for analyzing Rhodopsin genes in various species
           to build a relational tree based on Rhodopsin gene similarity.
           This ADT also builds a relational tree based on enviromental 
           variables for the same species, and then finally compares the two 
           trees to attempt to see if there is any relationship between 
           Rhodopsin gene similarites and environmental similarites.
           
Dependencies: BioPython

Assumptions and Implementation Notes:
            -- All dependencies must be accessible and configured correctly.
           

@authors: Camilo Acosta, Jayse Farrel, Jessica Kunder, Ryan Palm
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

__author__ = "Camilo Acosta, Jayse Farrel, Jessica Kunder, Ryan Palm"
__copyright__ = "COPYRIGHT_INFORMATION"
__credits__ = ["Camilo Acosta, Jayse Farrel, Jessica Kunder, Ryan Palm"]
__license__ = "GPL"
__version__ = "1.6.0dev"
__maintainer__ = "AUTHOR_NAME"
__email__ = "AUTHOR_EMAIL"
__status__ = "homework"

from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import Phylo
from Bio import AlignIO
from Bio.Alphabet import generic_dna
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor


__cephGenes = []
__sequences = []
__alignedGenes = []
__alignedGenesNorm = []
__numberOfGenes = 0 #starts at one for bootstrap case.

def main():

    #Parse input fasta file into list of strings containing genes
    infile = open ("sequence.txt", "r")
    parseInputFile(infile)
    
    print(__numberOfGenes)
    
    #Do alignments, grab alignments and splice to same length 
    #to generate alignment file
    processAlignments()
    
    #send alignment file to phyloTreeMaker to turn into a tree
    aln = AlignIO.read("asequences.txt", "fasta")
    phyloTreeMaker(aln)
    
    return

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
parseInputFile    : Reads input file and adds genes to global list
precondition      : input file must exist and be in FASTA format
postcondition     : genes added to __cephGenes
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''     
def parseInputFile(infile):
    
    global __cephGenes       # To modify global __avgGeneLength
    global __numberOfGenes
    
    tempGene = ''                # Holds the current gene
    
    tempLine = infile.readline() # boot strap case, remove first line
    __numberOfGenes +=1
    
    while True:
        
        tempLine = infile.readline()
        
        # Tail case (end of file reached)
        if tempLine == '':
            
            # final tempGene concatenate and add length to list then break
            strippedTempLine = tempLine.strip()
            tempGene += strippedTempLine
            
            #add gene to list here
            __cephGenes.append(tempGene)
            break
        
        if tempLine[0] == '>':
            
            __numberOfGenes += 1
            # On a heading line, tempGene done, add str and clear for next
            #add gene here list here
            __cephGenes.append(tempGene)
            
            tempGene = ''
        
        else:
            
            # building each gene line by line.
            strippedTempLine = tempLine.strip()
            tempGene += strippedTempLine


    
    return

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
processAlignments : align's genes and adds to matrix
precondition      : genes list must be initialized
postcondition     : matrix is updated with alignments scores
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''     
def processAlignments():
    #globals
    global __cephGenes
    global __numberOfGenes
    global __alignedGenes
    global __alignedGenesNorm
    global __numberOfGenes
    
    # run global alignments (right now only compairing the first gene to all 
    # others, will implement "to every other" functionality
    for x in range(0, __numberOfGenes):
     
        alignments = pairwise2.align.globalxx(__cephGenes[0], __cephGenes[x])
        
        __alignedGenes.append(alignments[0][1]) #add alignment to aligned genes
        
        print(pairwise2.format_alignment(*alignments[0]))
        print()
     
    
    print("Length of alignments before splicing:")
    for x in range(0, __numberOfGenes):
        print(len(__alignedGenes[x]))
    
    
    #Must splice genes based on the shortest alignment to create aligned file
    # (all genes must be the same length)
    for x in range(0,__numberOfGenes):
        temp = __alignedGenes[x]
        
        __alignedGenesNorm.append(temp[0:1360])
    
    # Need to redo below, hard way alignment creation
        
    
    #tempSeqList = []
    
    # bootstrap case, requre
    align1 = MultipleSeqAlignment([SeqRecord(Seq(__alignedGenesNorm[0], generic_dna,), id="Octopus bimaculoides")])
    
    # need to change "id=x" to hashmap of names
    for x in range(1, __numberOfGenes):
        align1.append(SeqRecord(Seq(__alignedGenesNorm[x], generic_dna,), id=''+str(x))) 
    
    
    #generate aligned file
    AlignIO.write(align1, "asequences.txt", "fasta")
    
    return


'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
phyloTreeMaker    : creates a phylogenetic tree from an aligned fasta file
precondition      : input file must be aligned
postcondition     : phylogenetic trees printed to output
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''      
def phyloTreeMaker(aln): 
    #make calculator object using default 'identiy' format (for dna)
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(aln)
    
    print(dm)
    print()
    
    #Create basic text tree and display
    constructor = DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.build_tree(aln)
    print(tree)
    print()
    
    # draw same tree using ascii art for tree structure
    Phylo.draw_ascii(tree)
    print()
    
    #generate image file (will pop up)
    Phylo.draw(tree)
    
    return
    
main()
