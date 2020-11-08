# Blast is lanched with the os module 
import os 

class BlastProcess: 

    # This counter is incremented when an instance is created 
    # basically we will include this counter in the naming of files in order to avoid overwritings
    counter = 0

    # constructor  
    def __init__(self, firstSequence, secondSequence, outformat):

        BlastProcess.counter += 1
        self.firstSequence = firstSequence
        self.secondSequence = secondSequence
        self.outformat = outformat 
        self.outfile = 'blast_' + str(BlastProcess.counter) + '_out_' + self.firstSequence[0:self.firstSequence.index('.')] + '_vs_' + self.secondSequence[0:self.secondSequence.index('.')] + '.txt'
        self.bestHitsFastaFile = 'best_hits_'+ str(BlastProcess.counter) + '_' + self.firstSequence[0:self.firstSequence.index('.')] + '_vs_' + self.secondSequence[0:self.secondSequence.index('.')] +'.fasta'
        self.reciprocalBlastFile = 'reciprocal_blast_'+ str(BlastProcess.counter) + '_' + self.firstSequence[0:self.firstSequence.index('.')] + '_vs_' + self.secondSequence[0:self.secondSequence.index('.')] + '.txt'
        # besthitsfastafile and reciprocalblastfile are written in this way to avoid overwriting

        # Create the file it doesn't exist
        file = open(self.outfile, 'w'); file.close()
    
    # This method is created for the blast process 
    def blastp(self):
        return os.system('blastp -query ' + self.firstSequence + ' -out ' + self.outfile + ' -subject ' + self.secondSequence + ' -outfmt ' + str(self.outformat) + ' -num_threads 2')        
        
    # This method will parse blast files 
    def parseBlastFile(self):
        # Import SearchIO module
        from Bio import SearchIO

        # Check the format of the format 
        if self.outformat == 5: 
            # Check if the outfile exists or not 
            if os.stat(self.outfile).st_size != 0: 
                # This for loop will parse in the generator and append it in the parsedFile list
                parsedFile = [i for i in SearchIO.parse(self.outfile, 'blast-xml')]

            # In case the outfile is empty we will proceed the blastprocess
            else: 
                self.blastp()
                # This for loop will parse in the generator and append it in the parsedFile list
                parsedFile = [i for i in SearchIO.parse(self.outfile, 'blast-xml')]
        
        # If the format is tabulated
        elif self.outformat == 6 or self.outformat == 7: 
            
            if os.stat(self.outfile).st_size != 0:
                # parsedFileGenerator is a generator it should be converted to a real list 
                parsedFileGenerator = SearchIO.parse(self.outfile, 'blast-tab', comments=True) if self.outformat == 7 else SearchIO.parse(self.outfile, 'blast-tab')
                # This for loop will parse in the generator and append it in the parsedFile list
                parsedFile = [i for i in parsedFileGenerator]

            else:
                self.blastp()
                parsedFileGenerator = SearchIO.parse(self.outfile, 'blast-tab', comments=True) if self.outformat == 7 else SearchIO.parse(self.outfile, 'blast-tab')
                # This for loop will parse in the generator and append it in the parsedFile list
                parsedFile = [i for i in parsedFileGenerator]

        # If the format is text 
        else: 
            if os.stat(self.outfile).st_size != 0:
                # This for loop will parse in the generator and append it in the parsedFile list
                parsedFile = [i for i in SearchIO.parse(self.outfile, 'blast-text')]

            else:
                self.blastp()
                # This for loop will parse in the generator and append it in the parsedFile list
                parsedFile = [i for i in SearchIO.parse(self.outfile, 'blast-text')]

        return parsedFile
    
    # This method is used to get the best hits 
    def getBestHits(self):

        parsedFile = self.parseBlastFile()
        # Open the fasta file which will contain the best hits 
        bestHitsFile = open(self.bestHitsFastaFile, 'w')

        for i in range(len(parsedFile)):
            # Create a list of HSPs
            HspFile = parsedFile[i].hsps

            # Now we will write the best hits in the fasta file 
            for j in range(len(HspFile)):
                if HspFile[j].evalue < 10 ** -100: 
                    # Write the best hits ids in the best_hits.fasta file
                    bestHitsFile.write("{}\n".format(HspFile[j].hit_id))

        bestHitsFile.close()

    # This method is used for the reciprocal blast
    def reciprocalBlast(self):
        # Check if the fasta file with the best hits exists or not 
        if os.path.isfile(self.bestHitsFastaFile):
            # launch the reciprocal blast 
            return os.system('blastp -query ' + self.bestHitsFastaFile + ' -out ' + self.reciprocalBlastFile + ' -subject ' + self.firstSequence + ' -outfmt 7 -max_target_seqs 1')
        
        else: 
            # otherwise if the fasta file doesn't exist, we will creat it first
            self.getBestHits()
            # then launch the reciprocal blast
            return os.system('blastp -query ' + self.bestHitsFastaFile + ' -out ' + self.reciprocalBlastFile  + ' -subject ' + self.firstSequence + ' -outfmt 7 -max_target_seqs 1')
              
# blastinstance = BlastProcess('Yersinia_pestis_angola.fasta', 'protéomes_yersia.fasta', 7)
 
# blastinstance.blastp()