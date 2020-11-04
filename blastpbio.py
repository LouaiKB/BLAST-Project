# Blast is lanched with the os module 
import os 

class BlastProcess: 

    # constructor  
    def __init__(self, firstSequence, secondSequence, outfile, outformat):
        # Create the file it doesn't exist
        file = open(outfile, 'w'); file.close()
        self.firstSequence = firstSequence
        self.secondSequence = secondSequence 
        self.outfile = outfile
        self.outformat = outformat

    # This method is created for the blast process 
    def blastp(self):
        return os.system('blastp -query ' + self.firstSequence + ' -out ' + self.outfile + ' -subject ' + self.secondSequence + ' -outfmt ' + str(self.outformat))        
        
    # This method will parse blast files 
    def parseBlastFile(self):
        # Import SearchIO module
        from Bio import SearchIO

        # Check the format of the format 
        if self.outformat == 5: 
            # Check if the outfile exists or not 
            if os.path.isfile(self.outfile): 
                # parsedFileGenerator is a generator it should be converted to a real list 
                parsedFileGenerator = SearchIO.parse(self.outfile, 'blast-xml')
                parsedFile = list()

                # This for loop will parse in the generator and append it in the parsedFile list
                for i in parsedFileGenerator:
                    # We will create a list with QueryResults
                    parsedFile.append(i)

            # In case the outfile is empty we will proceed the blastprocess
            else: 
                self.blastp()
                parsedFileGenerator = SearchIO.parse(self.outfile, 'blast-xml')
                parsedFile = list()

                # This for loop will parse in the generator and append it in the parsedFile list
                for i in parsedFileGenerator:
                    # We will create a list with QueryResults
                    parsedFile.append(i)
        
        # If the format is tabulated
        elif self.outformat == 6 or self.outformat == 7: 
            
            if os.path.isfile(self.outfile):
                # parsedFileGenerator is a generator it should be converted to a real list 
                parsedFileGenerator = SearchIO.parse(self.outfile, 'blast-tab', comments=True) if self.outformat == 7 else SearchIO.parse(self.outfile, 'blast-tab')

                parsedFile = list()

                # This for loop will parse in the generator and append it in the parsedFile list
                for i in parsedFileGenerator:
                    # We will create a list with QueryResults
                    parsedFile.append(i)                

            else:
                self.blastp()
                parsedFileGenerator = SearchIO.parse(self.outfile, 'blast-tab', comments=True) if self.outformat == 7 else SearchIO.parse(self.outfile, 'blast-tab')
                parsedFile = list()

                # This for loop will parse in the generator and append it in the parsedFile list
                for i in parsedFileGenerator:
                    # We will create a list with QueryResults
                    parsedFile.append(i)

        # If the format is text 
        else: 
            if os.path.isfile(self.outfile):
                # parsedFileGenerator is a generator it should be converted to a real list 
                parsedFileGenerator = SearchIO.parse(self.outfile, 'blast-text')
                parsedFile = list()

                # This for loop will parse in the generator and append it in the parsedFile list
                for i in parsedFileGenerator:
                    # We will create a list with QueryResults
                    parsedFile.append(i)                

            else:
                self.blastp()
                parsedFileGenerator = SearchIO.parse(self.outfile, 'blast-text')
                parsedFile = list()

                # This for loop will parse in the generator and append it in the parsedFile list
                for i in parsedFileGenerator:
                    # We will create a list with QueryResults
                    parsedFile.append(i)
        
        return parsedFile
    
    # This method is used to get the best hits 
    def getBestHits(self):

        parsedFile = self.parseBlastFile()
        # Open the fasta file which will contain the best hits 
        bestHitsFile = open('best_hits.fasta', 'w')

        for i in range(len(parsedFile)):
            # Create a list of HSPs
            HspFile = parsedFile[i].hsps

            # Now we will write the best hits in the fasta file 
            for j in range(len(HspFile)):
                if HspFile[j].evalue < 10 ** -50: 
                    if self.outformat == 5: 
                        # if the format is xml, we will write description and sequence
                        bestHitsFile.write(">{}\n".format(HspFile[j].hit_description))
                        bestHitsFile.write("{}\n".format(HspFile[j].hit.seq))
                    else:
                        # otherwise we will write only the id_hit
                        bestHitsFile.write("{}\n".format(HspFile[j].hit_id))
        bestHitsFile.close()

    # This method is used for the reciprocal blast
    def reciprocalBlast(self):
        # Check if the fasta file withe the best hits exists or not 
        if os.path.isfile('best_hits.fasta'):
            return os.system('blastp -query best_hits.fasta -out reciprocal_blast.txt -subject ' + self.secondSequence + ' -outfmt 7')
        
        else: 
            self.getBestHits()
            return os.system('blastp -query best_hits.fasta -out reciprocal_blast.txt -subject ' + self.secondSequence + ' -outfmt 7')
              