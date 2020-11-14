# import the BlastProcess Class 
from BlastpProcess import BlastProcess
from itertools import combinations 
from Bio import SearchIO
import pandas as pd 
import csv

# # store the files into variables
# yersiniaProteom = 'Yersinia_pestis_angola.fasta'
# yersiaProteom = 'protéomes_yersia.fasta'
# yersinaPestisBiovar = 'Yersinia_pestis_biovar_microtus_str_91001.ASM788v1.pep.all.fa'

# # lists of proteomes
# proteomeLists = [yersiniaProteom, yersiaProteom, yersinaPestisBiovar]

# # creation of combinations object 
# combinatioObject = combinations(proteomeLists, 2)

# # creation of proteom combinations 
# proteomeCombinations = [i for i in combinatioObject]

# # this list will store the files which have the results of reciprocal blast
# reciprocalBestHitsFiles = list()

# for i in range(len(proteomeCombinations)):
#     blastinstance = BlastProcess(proteomeCombinations[i][0], proteomeCombinations[i][1], 7)
#     blastinstance.getBestHits()
#     reciprocalBestHitsFiles.append(blastinstance.reciprocalBestHitsCsvFile)

reciprocalBestHitsFiles = ['reciprocal_best_hits/best_hits_2_Yersinia_pestis_angola_vs_protéomes_yersia.csv', 'reciprocal_best_hits/best_hits_3_Yersinia_pestis_angola_vs_Yersinia_pestis_biovar_microtus_str_91001.csv', 'reciprocal_best_hits/best_hits_4_protéomes_yersia_vs_Yersinia_pestis_biovar_microtus_str_91001.csv']
# Creation of variables 
all_results = dict()

# This for loop will retrieve ids of hits and querys from RBH files
for csv_file in reciprocalBestHitsFiles:

    # create a dateframe of the csv_files
    df = pd.read_csv(csv_file)
   
    # get the results (hit id, query id) from reciprocal files
    results = [i for i in df[['Hit_ID', 'Query_ID']].apply(tuple, axis=1)]
   
    # store the ids of hits and querys in a dictionnary
    all_results['results_'+ str(reciprocalBestHitsFiles.index(csv_file))] = results

# this function will return an Array with ids (cluster)
def clusterGenerator(all_results):
    cluster = dict()
   
    # retrieve values of dictionnary
    values = [i for i in all_results.values()]

    # this will loop will iterate throw values of the dictionnary (hits and querys)
    while len(values) > 0:
        concatenatedListOfValues = list()

        for j in range(1, len(values)):
            # create a concatenated list with the second to the last element of the values list
            concatenatedListOfValues += values[j]
        
        # this for loop will iterate throw elements of the first element of the values list 
        for result1 in values[0]:
            # this for loop will iterate throw elements of the concatenated list 
            for result2 in concatenatedListOfValues:
                # if there is a match between hits 
                if result1[0] == result2[0]:
                    # we will check if the result1[0] is a key in the cluster dictionnary 
                    if result1[0] in cluster.keys():
                        # this list contains the hits and the querys from the match 
                        resultsFromTheMatch = [result1[0], result1[1], result2[1]]
                        # we will parse the resultsFromTheMatch and check if these ids exist or not in the values of result1[0] key
                        for i in resultsFromTheMatch:
                            if (i) not in cluster[result1[0]]:
                                # if if doesn't exist we will concatenate it in the values which are tuples
                                cluster[result1[0]] += (i,)
                    
                    # if result1 doesn't exist in the keys of cluster dictionnary
                    else: 
                        resultsFromTheMatch = [result1[0], result1[1], result2[1]]
                        # this counter is used for checking
                        counter = 0 
                        # we will parse the resultsFromTheMatch 
                        for i in range(len(resultsFromTheMatch)):
                            # we will parse the keys in the cluster
                            for key in cluster.keys():
                                # we will check if the resultsFromTheMatch exists or not in the values to avoid redondance!
                                if (resultsFromTheMatch[i]) in cluster[key]:
                                    
                                    if i == 0:
                                        for z in (resultsFromTheMatch[1], resultsFromTheMatch[2]):
                                            # if the this ids don't exist in the values of the dictionnary add them 
                                            if z not in cluster[key]:
                                                cluster[key] += (z,); counter += 1

                                    elif i == 1:
                                        for z in (resultsFromTheMatch[0], resultsFromTheMatch[2]):
                                            if z not in cluster[key]:
                                                cluster[key] += (z,); counter += 1

                                    elif i == 2:
                                        for z in (resultsFromTheMatch[0], resultsFromTheMatch[1]):
                                            if z not in cluster[key]:
                                                cluster[key] += (z,); counter += 1
                        if counter == 0: 
                            # when the counter == 0 this means that all the resultsFromTheMatch don't exist in the dictionnary values
                            # that's why we create a new key in the cluster dictionnary and add the values in tuples
                            cluster[result1[0]] = (result1[0], result1[1], result2[1])
                
                # if there is a match between querys 
                if result1[1] == result2[1]:
                    # we will check if the result1[0] is a key in the cluster dictionnary 
                    if result1[1] in cluster.keys():
                        # this list contains the hits and the querys from the match 
                        resultsFromTheMatch = [result1[0], result1[1], result2[0]]
                        # we will parse the resultsFromTheMatch and check if these ids exist or not in the values of result1[0] key
                        for i in resultsFromTheMatch:
                            if (i) not in cluster[result1[1]]:
                                # if if doesn't exist we will concatenate it in the values which are tuples
                                cluster[result1[1]] += (i,)
                    
                    # if result1 doesn't exist in the keys of cluster dictionnary
                    else: 
                        resultsFromTheMatch = [result1[0], result1[1], result2[0]]
                        # this counter is used for checking
                        counter = 0 
                        # we will parse the resultsFromTheMatch 
                        for i in range(len(resultsFromTheMatch)):
                            # we will parse the keys in the cluster
                            for key in cluster.keys():
                                # we will check if the resultsFromTheMatch exists or not in the values to avoid redondance!
                                if (resultsFromTheMatch[i]) in cluster[key]:
                                    
                                    if i == 0:
                                        for z in (resultsFromTheMatch[1], resultsFromTheMatch[2]):
                                            # if the this ids don't exist in the values of the dictionnary add them 
                                            if z not in cluster[key]:
                                                cluster[key] += (z,); counter += 1

                                    elif i == 1:
                                        for z in (resultsFromTheMatch[0], resultsFromTheMatch[2]):
                                            if z not in cluster[key]:
                                                cluster[key] += (z,); counter += 1

                                    elif i == 2:
                                        for z in (resultsFromTheMatch[0], resultsFromTheMatch[1]):
                                            if z not in cluster[key]:
                                                cluster[key] += (z,); counter += 1
                        if counter == 0: 
                            # when the counter == 0 this means that all the resultsFromTheMatch don't exist in the dictionnary values
                            # that's why we create a new key in the cluster dictionnary and add the values in tuples
                            cluster[result1[1]] = (result1[0], result1[1], result2[0])
    
                # if there is a match between the hit in the first file and the query in the second file
                if result1[0] == result2[1]:
                    # we will check if the result1[0] is a key in the cluster dictionnary 
                    if result1[0] in cluster.keys():
                        # this list contains the hits and the querys from the match 
                        resultsFromTheMatch = [result1[0], result1[1], result2[0]]
                        # we will parse the resultsFromTheMatch and check if these ids exist or not in the values of result1[0] key
                        for i in resultsFromTheMatch:
                            if (i) not in cluster[result1[0]]:
                                # if if doesn't exist we will concatenate it in the values which are tuples
                                cluster[result1[0]] += (i,)
                    
                    # if result1 doesn't exist in the keys of cluster dictionnary
                    else: 
                        resultsFromTheMatch = [result1[0], result1[1], result2[0]]
                        # this counter is used for checking
                        counter = 0 
                        # we will parse the resultsFromTheMatch 
                        for i in range(len(resultsFromTheMatch)):
                            # we will parse the keys in the cluster
                            for key in cluster.keys():
                                # we will check if the resultsFromTheMatch exists or not in the values to avoid redondance!
                                if (resultsFromTheMatch[i]) in cluster[key]:
                                    
                                    if i == 0:
                                        for z in (resultsFromTheMatch[1], resultsFromTheMatch[2]):
                                            # if the this ids don't exist in the values of the dictionnary add them 
                                            if z not in cluster[key]:
                                                cluster[key] += (z,); counter += 1

                                    elif i == 1:
                                        for z in (resultsFromTheMatch[0], resultsFromTheMatch[2]):
                                            if z not in cluster[key]:
                                                cluster[key] += (z,); counter += 1

                                    elif i == 2:
                                        for z in (resultsFromTheMatch[0], resultsFromTheMatch[1]):
                                            if z not in cluster[key]:
                                                cluster[key] += (z,); counter += 1
                        if counter == 0: 
                            # when the counter == 0 this means that all the resultsFromTheMatch don't exist in the dictionnary values
                            # that's why we create a new key in the cluster dictionnary and add the values in tuples
                            cluster[result1[0]] = (result1[0], result1[1], result2[0])
                
                # if there is a match between query of the first file and hit of the second file 
                if result1[1] == result2[0]:
                    # we will check if the result1[0] is a key in the cluster dictionnary 
                    if result1[1] in cluster.keys():
                        # this list contains the hits and the querys from the match 
                        resultsFromTheMatch = [result1[0], result1[1], result2[1]]
                        # we will parse the resultsFromTheMatch and check if these ids exist or not in the values of result1[0] key
                        for i in resultsFromTheMatch:
                            if (i) not in cluster[result1[1]]:
                                # if if doesn't exist we will concatenate it in the values which are tuples
                                cluster[result1[1]] += (i,)
                    
                    # if result1 doesn't exist in the keys of cluster dictionnary
                    else: 
                        resultsFromTheMatch = [result1[0], result1[1], result2[1]]
                        # this counter is used for checking
                        counter = 0 
                        # we will parse the resultsFromTheMatch 
                        for i in range(len(resultsFromTheMatch)):
                            # we will parse the keys in the cluster
                            for key in cluster.keys():
                                # we will check if the resultsFromTheMatch exists or not in the values to avoid redondance!
                                if (resultsFromTheMatch[i]) in cluster[key]:
                                    
                                    if i == 0:
                                        for z in (resultsFromTheMatch[1], resultsFromTheMatch[2]):
                                            # if the this ids don't exist in the values of the dictionnary add them 
                                            if z not in cluster[key]:
                                                cluster[key] += (z,); counter += 1

                                    elif i == 1:
                                        for z in (resultsFromTheMatch[0], resultsFromTheMatch[2]):
                                            if z not in cluster[key]:
                                                cluster[key] += (z,); counter += 1

                                    elif i == 2:
                                        for z in (resultsFromTheMatch[0], resultsFromTheMatch[1]):
                                            if z not in cluster[key]:
                                                cluster[key] += (z,); counter += 1
                        if counter == 0: 
                            # when the counter == 0 this means that all the resultsFromTheMatch don't exist in the dictionnary values
                            # that's why we create a new key in the cluster dictionnary and add the values in tuples
                            cluster[result1[1]] = (result1[0], result1[1], result2[1])
        
        
        # now remove the first element in the value list and repeat the process 
        values.remove(values[0])
        
    # now write the cluster in a csv file
    clusterFile = open('cluster_all.csv', 'w')

    # initialise the writer of the csv 
    theWriter = csv.writer(clusterFile)

    # write in the file, we add list(set()) to avoid redondance
    for row in list(set(cluster.values())):
        theWriter.writerow(row)

    # close the file
    clusterFile.close()

clusterGenerator(all_results)

