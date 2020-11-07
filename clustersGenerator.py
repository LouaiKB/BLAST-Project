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
# reciprocalFiles = list()

# for i in range(len(proteomeCombinations)):
#     blastinstance = BlastProcess(proteomeCombinations[i][0], proteomeCombinations[i][1], 7)
#     blastinstance.reciprocalBlast()
#     reciprocalFiles.append(blastinstance.reciprocalBlastFile)

reciprocalFiles = [
    'reciprocal_blast_1_Yersinia_pestis_angola_vs_protéomes_yersia.txt', 
    'reciprocal_blast_2_Yersinia_pestis_angola_vs_Yersinia_pestis_biovar_microtus_str_91001.txt', 
    'reciprocal_blast_3_protéomes_yersia_vs_Yersinia_pestis_biovar_microtus_str_91001.txt'
]

all_results_ids = {}
all_results_querys = {}


for i in reciprocalFiles:
    queryResults = [i for i in SearchIO.parse(i, 'blast-tab', comments=True)]
    query_ids = [i[0].query_id for i in (i.hits for i in queryResults)]
    id_hits = [i[0].id for i in (i.hits for i in queryResults)]
    all_results_ids['hits_id_'+ str(reciprocalFiles.index(i)) + '_results'] = id_hits
    all_results_querys['query_id_' + str(reciprocalFiles.index(i)) + '_results'] = query_ids

print(all_results_ids)
print('-----------------------------------------')
print(all_results_querys)