
# launching blast using a system command
import os 

class BlastProcess: 

    # constructor
    def __init__(self, firstSequence, secondSequence, outfile, outformat):
        self.firstSequence = firstSequence
        self.secondSequence = secondSequence 
        self.outfile = outfile
        self.outformat = outformat

    # this method is used to do the blast process
    def blastp(self):
        return os.system('blastp -query ' + self.firstSequence + ' -out ' + self.outfile + ' -subject ' + self.secondSequence + ' -outfmt ' + str(self.outformat))

    # this method is used to get best hits for each protein 
    def getBestHits(self):
        # import of the module to parse the blast output file
        from Bio import SearchIO
        
        # outfmt 5 return a file in xml format 
        if self.outformat == 5: 
            
            #check if the file is empty or not
            if os.stat(self.outfile).st_size != 0:
                bestHitsFile = open('best_hits.fasta', 'w')
                # if the file is not empty, you don't need to blast it
                # parsedFileGenrator is a generator
                parsedFileGenerator = SearchIO.parse(self.outfile, 'blast-xml')
                parsedFile = list()

                # this for loop will parse in the generator created and do Append in the parsedFile list
                for i in parsedFileGenerator:
                    # create a list that contains the results of the QueryResult blast
                    parsedFile.append(i)

                # this for loop will get all the hits for each query
                for i in range(len(parsedFile)):
                    HspFile = parsedFile[i].hsps

                    # Write the best hits in the file according to the value condition of the e value
                    for j in range(len(HspFile)):
                        if HspFile[j].evalue < 10 ** -50: 
                            bestHitsFile.write(">{}\n".format(HspFile[j].hit_description))
                            bestHitsFile.write("{}\n".format(HspFile[j].hit.seq))

                bestHitsFile.close()
            else:
                # if the out file is empty, run the blast
                self.blastp()
                bestHitsFile = open('best_hits.fasta', 'w')
                parsedFileGenerator = SearchIO.parse(self.outfile, 'blast-xml')
                parsedFile = list()

                # this for loop will parse in the generator created and do Append in the parsedFile list
                for i in parsedFileGenerator:
                    # create a list that contains the results of the QueryResult blast
                    parsedFile.append(i)

                # this for loop will get all the hits for each query 
                for i in range(len(parsedFile)):
                    HspFile = parsedFile[i].hsps

                    # Write the best hits in the file according to the value condition of the e value
                    if HspFile.evalue < 10 ** -50: 
                        bestHitsFile.write(">{}\n".format(HspFile[j].hit_description))
                        bestHitsFile.write("{}\n".format(HspFile[j].hit.seq))
                
                bestHitsFile.close()


blast = BlastProcess('Yersinia_pestis_angola.fasta', 'proteìomes_yersia.fasta/protéomes_yersia.fasta', 'blast_out.xml', 5)

blast.getBestHits()
