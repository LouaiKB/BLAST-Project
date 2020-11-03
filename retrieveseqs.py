from Bio import Entrez, SeqIO

class ResearchFromNcbi: 
    
    saveNumb = 0 
    
    def __init__(self, email, research, database): 
        self.email = email 
        self.research = research 
        self.database = database  
        self.savingPath = "" 

    def doResearch(self): 
        try: 
            Entrez.email = self.email
            return Entrez.esearch(db=self.database, term=self.research, rettype='gb')
        except IOError:
            return "No connection!"

    def retrieve(self):
        try: 
            handle = self.doResearch()
            record = Entrez.read(handle)
            for i in range(len(record['IdList'])):
                newHandle = Entrez.efetch(db=self.database, id=record['IdList'][i], rettype='gb', retmode='text')
                newRecord = SeqIO.read(newHandle, 'gb')
                print("[",str(i + 1 ),"]", "Ident: " + newRecord.id + " -- Description: " + newRecord.description)
        except IOError: 
            return "No connection!"

    def saveOneFile(self):
        try:
            saveNumb = int(input('Which number: ')) 
            self.savingPath = input('Saving path: ')
            record = Entrez.read(self.doResearch())
            handle = Entrez.efetch(db=self.database, id=record['IdList'][saveNumb - 1], rettype='gb', retmode='text')
            record = SeqIO.read(handle, 'gb')
            SeqIO.write(record, self.savingPath + "\\" + record.id + ".gb", 'gb')
            print('successfully downloaded!')
        except IOError: 
            return "No connection!"

    def saveAllFiles(self):
        try:
            record = Entrez.read(self.doResearch())
            for i in range(len(record['IdList'])):
                handle = Entrez.efetch(db=self.database, id=record['IdList'][i], rettype='gb', retmode='text')
                record = SeqIO.read(handle, 'gb')
                SeqIO.write(record, self.savingPath + "\'" + record.id + ".gb", 'gb')
                print('successfully downloaded!')
        except IOError:
            return "No connection!"