# import the BlastProcess Class 
from blastpbio import BlastProcess
from itertools import combinations 
from Bio import SearchIO

# store the files into variables
yersiniaProteom = 'Yersinia_pestis_angola.fasta'
yersiaProteom = 'protéomes_yersia.fasta'
yersinaPestisBiovar = 'Yersinia_pestis_biovar_microtus_str_91001.ASM788v1.pep.all.fa'

# lists of proteomes
proteomeLists = [yersiniaProteom, yersiaProteom, yersinaPestisBiovar]

# creation of combinations object 
combinatioObject = combinations(proteomeLists, 2)

# creation of proteom combinations 
proteomeCombinations = [i for i in combinatioObject]

# this list will store the files which have the results of reciprocal blast
reciprocalFiles = list()

for i in range(len(proteomeCombinations)):
    blastinstance = BlastProcess(proteomeCombinations[i][0], proteomeCombinations[i][1], 7)
    blastinstance.reciprocalBlast()
    reciprocalFiles.append(blastinstance.reciprocalBlastFile)

# Creation of variables 
all_results_hits = {}
all_results_querys = {}
hits_cluster = list()

# This for loop will retrieve ids of hits and querys from RBH files
for i in reciprocalFiles:
    queryResults = [i for i in SearchIO.parse(i, 'blast-tab', comments=True)]
   
    # get the querys id from reciprocal files
    query_ids = [i[0].query_id for i in (i.hits for i in queryResults)]
   
    # get the hits id from recirocal files
    id_hits = [i[0].id for i in (i.hits for i in queryResults)]
   
    # store the ids of hits and querys in two dictionnaries
    all_results_hits['hits_id_'+ str(reciprocalFiles.index(i)) + '_results'] = id_hits
    all_results_querys['query_id_' + str(reciprocalFiles.index(i)) + '_results'] = query_ids

# this function will return an Array with ids (cluster)
def clusterGenerator(hits_id, querys_id):
    cluster = list()
   
    # retrieve values of dictionnaries in two variables
    all_hits = [i for i in hits_id.values()]
    all_querys = [i for i in querys_id.values()]

    # this will loop will iterate throw values of the dictionnarie (hits and ids)
    while len(all_hits) > 0 and len(all_querys) > 0:
        concatenatedListOfHits = list(); concatenatedListOfQuerys = list()

        for j in range(1, len(all_hits)):
            # create a concatenated list with the second to the last element of the all_hits list
            concatenatedListOfHits += all_hits[j]
        
        # this for loop will iterate throw elements of the first element of the all_hits list 
        for hitInTheFirstElement in all_hits[0]:
            # this for loop will iterate throw elements of the concatenated list 
            for hitInTheOtherElements in concatenatedListOfHits:
                # if there is a match between hits and the hits is not existed in the cluster then append the id in the cluster
                if hitInTheFirstElement == hitInTheOtherElements and hitInTheFirstElement not in cluster:
                    cluster.append(hitInTheFirstElement)
    
        for j in range(1, len((all_querys))):
            # create a concatenated list with the second to the last element of the all_querys list
            concatenatedListOfQuerys += all_querys[j]

        # this for loop will iterate throw elements of the first element of the all_querys list 
        for queryInTheFirstElement in all_querys[0]:
            # this for loop will iterate throw elements of the concatenated list 
            for queryInTheOtherElements in concatenatedListOfQuerys:
                # if there is a match between hits and the hits is not existed in the cluster then append the id in the cluster
                if queryInTheFirstElement == queryInTheOtherElements and queryInTheFirstElement not in cluster:
                    cluster.append(queryInTheFirstElement)

        # now remove the first element in these lists and repeat the process 
        all_hits.remove(all_hits[0])
        all_querys.remove(all_querys[0])
    
    return cluster
