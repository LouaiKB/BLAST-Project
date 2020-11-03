# On va lancer le blast en utilisant une commande system
import os 

class BlastProcess: 

    # constructeur 
    def __init__(self, firstSequence, secondSequence, outfile, outformat):
        self.firstSequence = firstSequence
        self.secondSequence = secondSequence 
        self.outfile = outfile
        self.outformat = outformat

    # cette méthode est utilisée pour faire le processus du blast 
    def blastp(self):
        return os.system('blastp -query ' + self.firstSequence + ' -out ' + self.outfile + ' -subject ' + self.secondSequence + ' -outfmt ' + str(self.outformat))

    # cette méthode est utilisée pour la récupération des meilleurs hits 
    def getBestHits(self):
        # importation du module qui permet de faire le parse de fichier de sortie du blast 
        from Bio import SearchIO
        
        # ce outfmt 5 return un fichier au format xml
        if self.outformat == 5: 
            
            # On va vérifier si le fichier est vide ou non 
            if os.stat(self.outfile).st_size != 0:
                bestHitsFile = open('best_hits.fasta', 'w')
                # Si le fichier n'est pas vide, c'est pas la peine de faire un blast
                # parsedFileGenrator est un générateur
                parsedFileGenerator = SearchIO.parse(self.outfile, 'blast-xml')
                parsedFile = list()

                # Cette boucle for va parser dans le générateur créer et faire Append dans la liste parsedFile
                for i in parsedFileGenerator:
                    # création d'une liste qui contient les résultats du blast QueryResult
                    parsedFile.append(i)

                # Cette boucle for va récupérer tous les hits pour chaque query 
                for i in range(len(parsedFile)):
                    HspFile = parsedFile[i].hsps

                    # Maintenant on va écrire les meilleurs hits qui ont un evalue 
                    for j in range(len(HspFile)):
                        if HspFile[j].evalue < 10 ** -50: 
                            bestHitsFile.write(">{}\n".format(HspFile[j].hit_description))
                            bestHitsFile.write("{}\n".format(HspFile[j].hit.seq))

                bestHitsFile.close()
            else:
                # Par contre si le fichier out file est vide il faut faire le blast  
                self.blastp()
                bestHitsFile = open('best_hits.fasta', 'w')
                parsedFileGenerator = SearchIO.parse(self.outfile, 'blast-xml')
                parsedFile = list()

                # Cette boucle for va parser dans le générateur créer et faire Append dans la liste parsedFile
                for i in parsedFileGenerator:
                    # création d'une liste qui contient les résultats du blast QueryResult
                    parsedFile.append(i)

                # Cette boucle for va récupérer tous les hits pour chaque query 
                for i in range(len(parsedFile)):
                    HspFile = parsedFile[i].hsps

                    # Maintenant on va écrire les meilleurs hits qui ont un evalue 
                    if HspFile.evalue < 10 ** -50: 
                        bestHitsFile.write(">{}\n".format(HspFile[j].hit_description))
                        bestHitsFile.write("{}\n".format(HspFile[j].hit.seq))
                
                bestHitsFile.close()


blast = BlastProcess('Yersinia_pestis_angola.fasta', 'proteìomes_yersia.fasta/protéomes_yersia.fasta', 'blast_out.xml', 5)

blast.getBestHits()