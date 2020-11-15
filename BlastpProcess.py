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
        self.outfile = 'blastp_out_files/blast_' + str(BlastProcess.counter) + '_out_' + self.firstSequence[0:self.firstSequence.index('.')] + '_vs_' + self.secondSequence[0:self.secondSequence.index('.')] + '.txt'
        self.outfile_reverse_blast = 'blastp_reverse_out_file/blast_reverse_' + str(BlastProcess.counter) + '_' + self.secondSequence[0:self.secondSequence.index('.')] + '_vs_' + self.firstSequence[0:self.firstSequence.index('.')] + '.txt'
        self.reciprocalBestHitsCsvFile = 'reciprocal_best_hits/best_hits_'+ str(BlastProcess.counter) + '_' + self.firstSequence[0:self.firstSequence.index('.')] + '_vs_' + self.secondSequence[0:self.secondSequence.index('.')] +'.csv'
        self.reciprocalBlastFile = 'reciprocal_files/reciprocal_blast_'+ str(BlastProcess.counter) + '_' + self.firstSequence[0:self.firstSequence.index('.')] + '_vs_' + self.secondSequence[0:self.secondSequence.index('.')] + '.txt'
        # reciprocalBestHitsCsvFile and reciprocalBlastFile are written in this way to avoid overwriting

        # Create the file it doesn't exist
        file = open(self.outfile, 'w'); file.close()
    
    # This method is created for the blast process 
    def blastp(self):
        return os.system('blastp -query ' + self.firstSequence + ' -out ' + self.outfile + ' -subject ' + self.secondSequence + ' -outfmt ' + str(self.outformat))        
        
    
    # This method is created for the reverse blast 
    def reverseBlastp(self):
        return os.system('blastp -query ' + self.secondSequence + ' -out ' + self.outfile_reverse_blast + ' -subject ' + self.firstSequence + ' -outfmt ' + str(self.outformat))
   
    # This method will parse blast files 
    def parseBlastFile(self):
        # Import SearchIO module
        from Bio import SearchIO

        # Check the format of the format 
        if self.outformat == 5: 
            # Check if the outfile exists or not 
            if os.stat(self.outfile).st_size != 0 and os.stat(self.outfile_reverse_blast).st_size != 0: 
                # This for loop will parse in the generator and append it in the parsedFile list and in the parsedReverseList
                parsedFile = [i for i in SearchIO.parse(self.outfile, 'blast-xml')]
                parsedReverseFile = [i for i in SearchIO.parse(self.outfile_reverse_blast, 'blast-xml')]

            # In case the outfile is empty we will proceed the blastprocess
            else: 
                self.blastp(); self.reverseBlastp()
                # This for loop will parse in the generator and append it in the parsedFile list and in the parsedReverseList
                parsedFile = [i for i in SearchIO.parse(self.outfile, 'blast-xml')]
                parsedReverseFile = [i for i in SearchIO.parse(self.outfile_reverse_blast, 'blast-xml')]

        # If the format is tabulated
        elif self.outformat == 6 or self.outformat == 7: 
            
            if os.stat(self.outfile).st_size != 0 and os.stat(self.outfile_reverse_blast).st_size != 0:
                # parsedFileGenerator and parsedFileReverseGenerator are generators which should be converted to a real lists 
                parsedFileGenerator = SearchIO.parse(self.outfile, 'blast-tab', comments=True) if self.outformat == 7 else SearchIO.parse(self.outfile, 'blast-tab')
                parsedFileReverseGenerator = SearchIO.parse(self.outfile_reverse_blast, 'blas-tab', comments=True) if self.outformat == 7 else SearchIO.parse(self.outfile_reverse_blast, 'blast-tab')
                # This for loop will parse in the generator and append it in the parsedFile list and in the parsedReverseList
                parsedFile = [i for i in parsedFileGenerator]
                parsedReverseFile = [i for i in parsedFileReverseGenerator]

            else:

                self.blastp(); self.reverseBlastp()
                parsedFileGenerator = SearchIO.parse(self.outfile, 'blast-tab', comments=True) if self.outformat == 7 else SearchIO.parse(self.outfile, 'blast-tab')
                parsedFileReverseGenerator = SearchIO.parse(self.outfile_reverse_blast, 'blast-tab', comments=True) if self.outformat == 7 else SearchIO.parse(self.outfile_reverse_blast, 'blast-tab')
                # This for loop will parse in the generator and append it in the parsedFile list and in the parsedReverseList
                parsedFile = [i for i in parsedFileGenerator]
                parsedReverseFile = [i for i in parsedFileReverseGenerator]

        # If the format is text 
        else: 
            
            if os.stat(self.outfile).st_size != 0 or os.stat(self.outfile_reverse_blast).st_size != 0:
                # This for loop will parse in the generator and append it in the parsedFile list and in the parsedReverseList
                parsedFile = [i for i in SearchIO.parse(self.outfile, 'blast-text')]
                parsedReverseFile = [i for i in SearchIO.parse(self.outfile_reverse_blast, 'blast-text')]
            
            else:
                self.blastp(); self.reverseBlastp()
                # This for loop will parse in the generator and append it in the parsedFile list and in the parsedReverseList
                parsedFile = [i for i in SearchIO.parse(self.outfile, 'blast-text')]
                parsedReverseFile = [i for i in SearchIO.parse(self.outfile_reverse_blast, 'blast-text')]

        # these two lists are created for storing the results of the best hits
        resultsBlast = []; resultReverseBlast = []
        
        # this for loop is to parse in the queryResults of the out blast file in order to get the best hits
        for i in range(len(parsedFile)):
            # get the best hits from each query result
            besthit = parsedFile[i].hits
            
            # this try except block is used to prevent the IndexErrors 
            try:
                # we get the hit id and the query id of the best hit in the first blast
                resultsBlast.append((besthit[0].id, besthit[0].query_id))
            
            except IndexError: 
                # in case of IndexError continue the parsing and generating of best hits
                pass

        # this for loop is to parse in the queryResults of the out blast reverse file in order to get the best hits
        for i in range(len(parsedReverseFile)):
            # get the best hits from each query result
            besthit = parsedReverseFile[i].hits

            # this try except block is used to prevent the IndexErrors
            try: 
                # we get the query id and the hit id of the best hits in the reverse blast, note that we inversed the ordre of the tuple in order to 
                # make the comparaison easier
                resultReverseBlast.append((besthit[0].query_id, besthit[0].id))

            except IndexError: 
                 # in case of IndexError continue the parsing and generating of best hits
                pass
        
        # this method return a tuple of the best hits in the first and the second blast
        return (resultsBlast, resultReverseBlast)
    
    # This method is used to get the best hits 
    def getBestHits(self):
        # get the best hits of the two blasts (tuple)
        parsedFile = self.parseBlastFile()
        
        # get the results of tuple in two variables which are lists now         
        bestHitsOfTheFirstBlast = parsedFile[0]
        bestHitsOfTheReverseBlast = parsedFile[1]

        # we create the output file which will contain the RBH
        rbhFile = open(self.reciprocalBestHitsCsvFile, 'w', newline='')
        
        # import the csv module  
        import csv
        # create an csv.DictWriter Object
        theWriter = csv.DictWriter(rbhFile, fieldnames=['Hit_ID', 'Query_ID'])
        
        # initialise the header
        theWriter.writeheader()

        # this for loop is used to parse the besthits of the first blast
        for itemInTheFirstBlast in bestHitsOfTheFirstBlast:
            # this for loop is used to parse the besthits of the reverse blast
            for itemInTheSecondBlast in bestHitsOfTheReverseBlast:
                    if itemInTheFirstBlast == itemInTheSecondBlast:
                        theWriter.writerow({
                            'Hit_ID': itemInTheFirstBlast[0],
                            'Query_ID': itemInTheFirstBlast[1]
                        })
        rbhFile.close()

    # This method is used for the blast with the best hits
    def blastWithBestHits(self):
        # Check if the fasta file with the best hits exists or not 
        if os.path.isfile(self.reciprocalBestHitsCsvFile):
            try:
                # launch the reciprocal blast 
                return os.system('blastp -query ' + self.reciprocalBestHitsCsvFile + ' -out ' + self.reciprocalBlastFile + ' -subject ' + self.firstSequence + ' -outfmt 7 -max_target_seqs 1')
            
            except RuntimeError:
                print('reciprocal blast done!')
        else: 
            # otherwise if the fasta file doesn't exist, we will creat it first
            self.getBestHits()
            # then launch the reciprocal blast
            return os.system('blastp -query ' + self.reciprocalBestHitsCsvFile + ' -out ' + self.reciprocalBlastFile  + ' -subject ' + self.firstSequence + ' -outfmt 7 -max_target_seqs 1')
              