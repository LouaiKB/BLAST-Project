from Bio import SeqIO, Entrez, Phylo
import pandas as pd 
import csv
import os
# store the files into variables
yersiniaProteom = 'Fervidicoccus.fasta'
yersiaProteom = 'Aciduliprofundum.fasta'
yersinaPestisBiovar = 'Ignicoccus.fasta'

# creation of a proteom list
proteomList = [yersiniaProteom, yersiaProteom, yersinaPestisBiovar]

ids_dict = dict()

# get the ids of proteins in a dictionnary ids_dict
for protein in proteomList:
    # the ids of each protein is stored in this dictionnary
    ids_dict['ids_of_' + protein[0:protein.index('.')]] = [i.id for i in (i for i in SeqIO.parse(protein, 'fasta'))]

# this function is used for sorting cluster
def sortClusters():
    # ids_dict is a global variable
    global ids_dict
    # create a list of all clusters
    clusters = pd.read_csv('cluster_all.csv', error_bad_lines=False).values.tolist()
    
    # this for loop is used to iterate throw the list of clusters
    for cluster in clusters:
        # this for loop is used to iterate throw the keys of the ids dictionnary 
        for key in ids_dict:
            counter = 0
            # this for loop is used to iterate throw each list of cluster
            for item in cluster: 
                # check if this item exists or not in the proteom 
                if item in ids_dict[key]:
                    counter += 1
            # if the counter > 1 this means that the same protein exists in two proteoms
            if counter > 1:
                # in this case we remove the cluster
                clusters.remove(cluster)

    # create a csv file which will contain the sorted clusters
    sortedFile = open('cluster_all_sorted.csv', 'w')
    # initialise the writer of the csv 
    theWriter = csv.writer(sortedFile)
    # we write each cluster list in the csv file
    for row in clusters:
        theWriter.writerow(row)
    # close the sortedFile
    sortedFile.close()

def alignementProcessOfClusters():
    global proteomList
    all_proteins = []
    for protein in proteomList:
        record = [i for i in SeqIO.parse(protein, 'fasta')]
        for i in record:
            all_proteins.append((i.id, i.description, i.seq))
    
    clusters = pd.read_csv('cluster_all_sorted.csv').values.tolist()
    for cluster in clusters:
        fastafile = open('infile.fasta', 'w')
        for item in cluster:
            for seqiorecord in all_proteins:
                if item in seqiorecord:
                    fastafile.write('>{}\n'.format(seqiorecord[0] + ' ' + seqiorecord[1]))
                    fastafile.write('{}\n'.format(seqiorecord[2]))
        fastafile.close()
        
        os.system('muscle -in infile.fasta -out out/outfile' + str(clusters.index(cluster)) +'.afa')

def createSuperAlignementsFile():
    from glob import glob
    from collections import defaultdict

    global proteomList
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
            for record in idAndSequence:
                if record[0] in all_proteins_ids[key]:
                    newProteinSeqs[key].append(record[1])
    
    superAlignementFile = open('superAlignFile.fasta', 'w')

    for key in newProteinSeqs:
        superAlignementFile.write('>{}\n'.format(key))
        for seq in newProteinSeqs[key]:
            superAlignementFile.write('{}'.format(seq))
        superAlignementFile.write('\n\n')
    
    superAlignementFile.close()

print(createSuperAlignementsFile())