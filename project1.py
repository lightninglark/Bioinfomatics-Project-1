'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
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

def main():

	human = open("FASTAs/Others/human.fasta", "r").read()
	redjunglefowl =  open("FASTAs/Others/redjunglefowl.fasta", "r").read()

	alignments = pairwise2.align.globalxx(human, redjunglefowl)
	print(alignments)