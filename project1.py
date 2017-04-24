'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
from tables.idxutils import ccs_ultralight
from tables.idxutils import ccs_ultralight
from Bio.FSSP.fssp_rec import align
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


__cephGenes = []
__alignmentMatrix = [[0]*3 for i in range(3)] #magic number 4, is the number of files we're dealing with

def main():
    
    infile = open ("sequence.txt", "r")
    parseInputFile(infile)
    #print(__cephGenes)
    
    processAlignments()
    

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
parseInputFile    : Reads input file and adds genes to global list
precondition      : input file must exist and be in FASTA format
postcondition     : genes added to __cephGenes
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''     
def parseInputFile(infile):
    
    global __cephGenes       # To modify global __avgGeneLength
    
    tempGene = ''                # Holds the current gene
    
    tempLine = infile.readline() # boot strap case, remove first line
    
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
    global __alignmentMatrix
    
    
    #test alignment
    alignments = pairwise2.align.globalxx(__cephGenes[0],__cephGenes[1])
    
    print(alignments[0])
    
    return
    
    
    
main()