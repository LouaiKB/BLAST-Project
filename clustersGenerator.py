# import the BlastProcess Class
from blastpbio import BlastProcess
from itertools import combinations 

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

for i in range(len(proteomeCombinations)):
    blastinstance = BlastProcess(proteomeCombinations[i][0], proteomeCombinations[i][1], 7)
    blastinstance.reciprocalBlast()