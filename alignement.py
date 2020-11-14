from Bio import SeqIO
import pandas as pd 
import csv

# store the files into variables
yersiniaProteom = 'Yersinia_pestis_angola.fasta'
yersiaProteom = 'protéomes_yersia.fasta'
yersinaPestisBiovar = 'Yersinia_pestis_biovar_microtus_str_91001.ASM788v1.pep.all.fa'

# creation of a proteom list
proteomList = [yersiniaProteom, yersiaProteom, yersinaPestisBiovar]

ids_dict = dict()

# get the ids of proteins in a dictionnary ids_dict
for protein in proteomList:
    # the ids of each protein is stored in this dictionnary
    ids_dict['ids_of_' + protein[0:protein.index('.')]] = [i.id for i in (i for i in SeqIO.parse(protein, 'fasta'))]


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
    
    sortedFile.close()

sortClusters()