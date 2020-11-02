# Importation de librairie BioPython 
# NcbiblastnCommandline est une classe pour faire Nucléotide-Nucléotide BLAST 
# NcbiblastpCommandline est une classe pour faire Protein-Protein BLAST
# More explications on https://biopython.org/docs/1.75/api/Bio.Blast.Applications.html

from Bio.Blast.Applications import NcbiblastnCommandline, NcbiblastpCommandline 

# nucleotideBlast est une fonction pour faire le blast nucléotide-nucléotide
def nucleotideBlast(query, cible, outfile):
    # Création d'un objet de la class NcbiblastnCommandline
    # outfmt 6 pour le format tabulé, outfmt 5 pour le format xml 
    blastProcess = NcbiblastnCommandline(query=query, subject=cible, outfmt=6, out=outfile)



