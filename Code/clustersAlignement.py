"""

AUTHOR : Louai KASSA BAGHDOUCHE

"""
from Bio import SeqIO, Entrez, Phylo, AlignIO
from glob import glob
from collections import defaultdict
import os


# creation of a proteom list
proteoms = input('Enter the proteoms with this format (proteom1,protem2, ...): ')

proteomList = proteoms.split(',')

ids_dict = dict()

# get the ids of proteins in a dictionnary ids_dict
for protein in proteomList:
    # the ids of each protein is stored in this dictionnary
    ids_dict['ids_of_' + protein[0:protein.index('.')]] = [i.id for i in (j for j in SeqIO.parse(protein, 'fasta'))]

# this function is used for sorting cluster
def sortClusters():
    # ids_dict is a global variable
    global ids_dict, proteomList
    # create a list of all clusters
    with open('cluster_all.txt', 'r') as f:
        lines = f.readlines()

    # this new file will have the sorted clusters
    with open('cluster_all_sorted.txt', 'w') as file:
        for line in lines:
            # create list of each cluster
            cluster = line.rstrip('\n').split(',')
            # 
            if len(cluster) <= len(proteomList):
                counter = 0 
                for key in ids_dict:
                    for item in cluster:
                        if item in ids_dict[key]:
                            counter += 1
                if counter <= len(ids_dict): 
                    file.write(','.join(cluster)+"\n")

# this function will proceed alignements between files 
def alignementProcessOfClusters():
    
    global proteomList
    # this for loop will create a list with all the ids and description and the sequences of the proteoms
    all_proteins = []
    for protein in proteomList:
        record = [i for i in SeqIO.parse(protein, 'fasta')]
        for i in record:
            all_proteins.append((i.id, i.description, i.seq))
    
    # read the clusters files (sorted)
    with open('cluster_all_sorted.txt', 'r') as f:
        lines = f.readlines()
    
    # this for loop is created for the alignments of clusters 
    for line in lines:
        cluster = line.rstrip('\n').split(',')
        fastafile = open('infile.fasta', 'w')
        for item in cluster:
            for seqiorecord in all_proteins:
                if item in seqiorecord:
                    fastafile.write('>{}\n'.format(seqiorecord[0] + ' ' + seqiorecord[1]))
                    fastafile.write('{}\n'.format(seqiorecord[2]))
        fastafile.close()
        
        os.system('muscle -in infile.fasta -out out/outfile' + str(lines.index(line)) +'.afa')

# function to generate the super Alignment file
def createSuperAlignementsFile():

    global proteomList

    # this for loop will generate a dictionnary of proteoms with there ids
    all_proteins_ids = {}
    for protein in proteomList: 
        record = [i for i in SeqIO.parse(protein, 'fasta')]
        all_proteins_ids[protein[0:protein.index('.')]] = []
        for i in record:
            all_proteins_ids[protein[0:protein.index('.')]].append(i.id)
    
    newProteinSeqs = defaultdict(list)
    
    for afaFile in glob('out/*.afa'):
        idAndSequence = [(i.id, i.seq) for i in SeqIO.parse(afaFile, 'fasta')]
        for key in all_proteins_ids.keys():
            if len(idAndSequence) == len(all_proteins_ids):
                for record in idAndSequence:
                    if record[0] in all_proteins_ids[key]:
                        newProteinSeqs[key].append(record[1])
            else:
                keys = [i for i in all_proteins_ids.keys()]
                for record in idAndSequence:
                    if record[0] in all_proteins_ids[key]:
                        newProteinSeqs[key].append(record[1])
                        keys.remove(key)
                    
                for key in keys:
                    newProteinSeqs[key].append('-' * len(idAndSequence[0][1]))


    superAlignementFile = open('superAlignFile.fasta', 'w')

    for key in newProteinSeqs:
        superAlignementFile.write('>{}\n'.format(key))
        for seq in newProteinSeqs[key]:
            superAlignementFile.write('{}'.format(seq))
        superAlignementFile.write('\n\n')
    
    superAlignementFile.close()

def generatePhylogeneticTree():
   # converting super alignment fasta file to phylip file
   fastaFile = AlignIO.parse('superAlignFile.fasta', 'fasta')
   phylipFile = AlignIO.write(fastaFile, 'superAlign.phy', 'phylip')
   
   # phylogenetic tree generator 
   tree = Phylo.read(phylipFile, 'newick')
   Phylo.draw(tree)
