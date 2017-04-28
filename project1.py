'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
from Crypto.PublicKey.DSA import construct
from Tests.test_TreeConstruction import NNITreeSearcherTest
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
            -- input file "sequence.txt" must be in FASTA format


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
from Bio.Phylo.TreeConstruction import DistanceCalculator, ParsimonyScorer,\
    NNITreeSearcher, ParsimonyTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.SubsMat.MatrixInfo import blosum62
from Bio.Phylo.Consensus import *


__sequences = []        # holds raw sequences from file
__alignedGenes = []     # holds aligned sequences
__alignedGenesNorm = [] # holds spliced (snipped) aligned gene sequences
__speciesNames = []     # holds the name of each specie
__numberOfGenes = 0     # holds number of genes

def main():

    #Parse input fasta file into list of strings containing genes
    infile = open ("sequence.txt", "r")
    parseInputFile(infile)

    # Do alignments, grab alignments and splice to same length
    # to generate alignment file
    processAlignments()

    #send alignment file to phyloTreeMaker to turn into a tree
    aln = AlignIO.read("asequences.txt", "fasta")
    phyloTreeMaker(aln)

    return

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
parseInputFile    : Reads input file and adds genes to global list
precondition      : input file must exist and be in FASTA format
postcondition     : genes added to __sequences, __numberOfGenes incremented
                    by one for each gene added, species name added to
                    __speciesNames
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
def parseInputFile(infile):

    # Forward reference to globals
    global __sequences
    global __numberOfGenes

    tempGene = ''                # Holds the current gene

    # Boot strap case (first gene case)
    tempLine = infile.readline() # grab first line
    speciesNameParser(tempLine)  # Parse species name and add to global list
    __numberOfGenes +=1          # add one gene/species

    while True:

        tempLine = infile.readline()

        # Tail case (end of file reached)
        if tempLine == '':

            # final tempGene concatenate and add length to list then break
            strippedTempLine = tempLine.strip()
            tempGene += strippedTempLine

            #add gene to list here
            __sequences.append(tempGene)
            break

        if tempLine[0] == '>':
            speciesNameParser(tempLine)  # Parse species name
            __numberOfGenes += 1
            # On a heading line, tempGene done, add str and clear for next
            #add gene here list here
            __sequences.append(tempGene)
            tempGene = ''

        else:
            # building each gene line by line.
            strippedTempLine = tempLine.strip()
            tempGene += strippedTempLine

    return

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
processAlignments : parses species name from fasta header line
precondition      : headerline must be formated like ">Species_Name|..."
postcondition     : species name is added to __speciesName list
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
def speciesNameParser(headerLine):

    # Forward reference
    global __speciesNames

    counter = 0 #counts length of name (in characters)

    #everything after the < token until the first pipe
    tempSpecieName = headerLine[1:headerLine.find('|')]

    __speciesNames.append(tempSpecieName)
    return

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
processAlignments : does gene alignment (comparison to a reference (first
                  : gene in inputfile (sequences.txt)))
precondition      : no significant precondition
postcondition     : aligned and "trimmed" genes are added to asequences.txt
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
def processAlignments():

    # Forward reference
    global __sequences
    global __numberOfGenes
    global __alignedGenes
    global __alignedGenesNorm
    global __numberOfGenes
    global __speciesNames

    # Run global alignments (first gene acts as a reference gene)
    for x in range(0, __numberOfGenes):

        alignments = pairwise2.align.globalds(__sequences[0], __sequences[x], blosum62, -10., -0.5)

        __alignedGenes.append(alignments[0][1]) #add alignment to aligned genes

        print(__speciesNames[x] + ": ")
        print(pairwise2.format_alignment(*alignments[0]))
        print()

    # Find shortest length of all aligned genes (all genes must be same
    # character length)
    shortestLength = len(min(__alignedGenes, key=len))
    
    print("Shortest Length Gene:" + str(shortestLength))

    # trim all genes based on shortest length variable
    for x in range(0,__numberOfGenes):
        temp = __alignedGenes[x]

        __alignedGenesNorm.append(temp[0:shortestLength])

    # Create MultipleSeqAlignment needed for creation of distance map
    # (and phylogenetic tree)

    # Bootstrap case (required by MultipleSeqAlignment constructor)
    align1 = MultipleSeqAlignment([SeqRecord(Seq(__alignedGenesNorm[0],
                                   generic_dna,), id=__speciesNames[0])
                                   ])

    # Append each additional alignment to the MSA
    for x in range(1, __numberOfGenes):
        align1.append(SeqRecord(Seq(__alignedGenesNorm[x], generic_dna,),
                                id=__speciesNames[x]))
    #generate aligned file
    AlignIO.write(align1, "asequences.txt", "fasta")

    return

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
phyloTreeMaker    : creates a phylogenetic tree from an aligned fasta file
precondition      : input file must exist and be aligned (and trimmed)
postcondition     : Distance map and Phylogenetic Tree structures
                    printed to output
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
    # Distance tree
    Phylo.draw(tree, show_confidence=True)


    #Parsimony Tree
    scorer = ParsimonyScorer()
    searcher = NNITreeSearcher(scorer)

    constructorParse = ParsimonyTreeConstructor(searcher, tree)
    parse_tree = constructorParse.build_tree(aln)
    Phylo.draw(parse_tree, show_confidence=True)

    '''
    # This tree works, but is EXTREMELY slow. I'd suggest adjusting the
    # lower if it is taking too long. 100 seems to take 5-10 minutes,
    # but the higher the replication number, the better the results.

    #Bootstrap (consensus trees)

    replicationNumber = 100

    calculatorCon = DistanceCalculator('blosum62')
    constructorCon = DistanceTreeConstructor(calculatorCon)
    treeCon = bootstrap_trees(aln, replicationNumber, constructorCon)

    majority_tree = majority_consensus(treeCon, 0.5)
    Phylo.draw(majority_tree, show_confidence=True)
    '''

    return

main()
