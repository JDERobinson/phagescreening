import re
import pandas as pd
from Bio import SeqIO

class processor(): #this processes a single .gb and fasta file with one entry 

            
    def __init__(self, fastaSeqs, gbFiles = None):
        self.genomes = self.__process_fasta(fastaSeqs) #a dictionary of all the genomes in a fasta file
       
        if gbFiles:    
           # print('I HAVE FOUND A GBFILE')
            self.proteomes = self.__process_gb(gbFiles) 
        else:
            self.proteomes = None
    
    ##methods##
    
    def __process_fasta(self, fastaSeqs):
        seqs = {}
        for record in SeqIO.parse(fastaSeqs,'fasta'):
            seqs[record.id] = str(record.seq)
        return seqs

    
    def __process_gb(self, gbFiles):
        
        prots = {}
   
        for record in SeqIO.parse(gbFiles, 'gb'):
                        
            feats = record.features
            protID = record.id
            
            PDTS = []
            TAGS = []
            LOCS = [] 
            
            for i in range(0, len(feats), 2):
                pdt = feats[i].qualifiers.get('product')
                PDTS.append(pdt) 

            for i in range(1, len(feats), 2):
                tag = feats[i].qualifiers['locus_tag'][0]
                TAGS.append(tag)
                loc = (int(feats[i].location.start) - 1, int(feats[i].location.end) -1) #get start and end positions within fasta, 
                                                                                        #python is 0 based indexing so gb start/end need to subtract 1
                LOCS.append(loc)
            
            INFO = list(zip(PDTS[1:], LOCS))               
            cdsRegions = dict(zip(TAGS, INFO))
            prots[protID] = cdsRegions
        
        return prots
    
        
    def __str__(self):
        if self.proteomes == None:
            genoPrint = 'Printing genomes...\n'
            for genome in self.genomes:
                genoPrint = f'Genome ID: {genome}\nGenome: {self.genomes[genome]}'
            return genoPrint
            
        else:
            genoPrint = 'Printing genomes...\n'
            for genome in self.genomes:
                genoPrint = f'Genome ID: {genome}\nGenome: {self.genomes[genome]}'
            
            protPrint = 'Proteome ID: '
            for proteomeID in self.proteomes:
                protPrint = protPrint + proteomeID + '\nProteins: '
                for protein in self.proteomes[proteomeID]:
                   # print(self.proteomes)
                    proteinLoc = self.proteomes[proteomeID][protein][0]
                    proteinPdt = self.proteomes[proteomeID][protein][1]
                    protPrint = protPrint + f'{protein}: ({proteinLoc.start}, {proteinLoc.end}) \nproduct: {proteinPdt}\n' 
            
            return f'{genoPrint}\n\n\n\n{protPrint}'        

class compiler(): #chunks fasta sequences by coding regions and returns organism objects
        
    def __init__(self, genomes, proteomes = None): #where proteomes and genomes are dictionaries
        self.__ids = genomes.keys()
        self.__genomes = genomes
        if proteomes:
            self.__proteomes = proteomes
            self.__organisms = self.__compile_organisms(genomes = self.__genomes, proteomes = self.__proteomes)
        else:
            self.__organisms = self.__compile_organisms(genomes = self.__genomes)
    
    
    class organism(): #match self.genomic to self.proteomic based on id
       
        def __init__(self, ID, orgGenome, orgProteome = None):
            self.id = ID
            self.genome = orgGenome
            
            if orgProteome:
                self.proteome = self.__chunk_CDS(self.genome, orgProteome)
                self.nonCDS = self.__chunk_NonCDS(self.genome, orgProteome)
            
            else: 
                self.proteome = 'Proteome not supplied'
                self.nonCDS = 'Proteome not supplied'

            
        
        def __chunk_CDS(self, cfGenome, cfProteome):
            CDSchunks = {}

            for protein in cfProteome:
               # print(cfProteome[protein])
                cds = []
                start = cfProteome[protein][1][0]
                end = cfProteome[protein][1][1] 
                pdt = cfProteome[protein][0]
              #  print(pdt)
                cds = [pdt, (start+1, end+1), cfGenome[start:end]] #python is zero based indexing, so had to substract 1 to slice correctly, 
                                                                    #location within genome is not actually zero based, so need to add back in 1 to get the 
                                                                    # right location
                CDSchunks[protein] = cds 
            return CDSchunks
        
        def __chunk_NonCDS(self, cfGenome, cfProteome): #chunk non-coding regions
            nonCDSChunks = {}
            TAG_LIST = [key for key in cfProteome.keys()]
            LENGTH = len(TAG_LIST)
        #    print(LENGTH)
            
            for i in range(LENGTH-1):

                if i == 0: #if at start of protein list
                    first = cfProteome[TAG_LIST[0]] #get first protein
                    firstStart = first[1][0] #get start location of first protein
                    nonCDSChunks[(0, firstStart+1)] = cfGenome[0:firstStart]
              #      print('this is the start of the list', i, firstStart)

                    nonCDSChunks = self.__help__chunk_NonCDS(i, TAG_LIST, cfProteome, cfGenome, nonCDSChunks)
                
                elif i == LENGTH-1: #if at the last gene of genome
                    last = cfProteome[TAG_LIST[LENGTH-1]]
                    lastEnd = last[1][1]
             #       print('this is the end of the list', i, lastEnd)
                    nonCDSChunks[lastEnd+1,len(cfGenome)] = cfGenome[lastEnd:len(cfGenome)-1]
                
                elif i < LENGTH-1:
                    nonCDSChunks = self.__help__chunk_NonCDS(i, TAG_LIST, cfProteome, cfGenome, nonCDSChunks)
        
            return nonCDSChunks


        def __help__chunk_NonCDS(self, i, TAG_LIST, cfProteome, cfGenome, nonCDSChunks):
            curr = cfProteome[TAG_LIST[i]] 
           # print('this is the current protein', curr)
            next = cfProteome[TAG_LIST[i+1]]
          #  print('this is the next protein', i+1, next)

            currEnd = curr[1][1] #get end location of current protein
            nextStart = next[1][0] #get start location of next protein

            if currEnd < nextStart:
                nonCDSChunks[(currEnd+1, nextStart+1)] = cfGenome[currEnd:nextStart] #get region of genome between end of curr to begininning of next
                                                                                    #stored by actual location in genome
            return nonCDSChunks

                        
        def __str__(self): ##summary statement of an organism, prints ID, length of genome, and number of genes + tags
            if self.proteome:
                protein_products = [self.proteome[item][0][0] for item in self.proteome.keys()]
                s = (f'Reference ID: {self.id}\nGenome length: {len(self.genome)}\nNumber of proteins: {len(self.proteome)}' 
                    f'\nProteins: {protein_products}')           
                return s
            
   
    def __compile_organisms(self, genomes, proteomes): #makes dictionary of all organisms
       
        organisms = {}

        for ID in genomes.keys():
            current_organism = self.organism(ID, genomes.get(ID), proteomes.get(ID))
            organisms[current_organism.id] = current_organism
    
        return organisms
  
        
    def get_organism(self, organismID): #allows user to access an organism object based on ID
        return self.__organisms[organismID]

    def all_organisms(self): #returns all organism objects as dictionaries, sorted by ID
        return self.__organisms
    

    ##TO DO
    #write function/class that can take an organism object and save the data to a text or excel file with protein:sequence 


class GPSgRNA():
    
    def __find_sgRNA(self): #helper that finds sgRNAs in a single sequence
        #define pams here 
        pass
    
    def get_chunked_sgRNA(self): #get all sgRNAs from a chunked fasta file -- calls find_sgRNA to find it
        pass
    
    
    def get_fasta_sgRNA(self): #get sgRNA's in fasta sequence
        pass


if __name__ == '__main__': 
    p = processor(fastaSeqs = 'Kampy.fasta', gbFiles = 'Kampy_annotated.gb')
   # print(p)

    c = compiler(p.genomes, p.proteomes)
   # print(c.organisms['KJ510414.1'].chunkedProteome)
   # print(c.all_organisms())
   # print(c.get_organism('KJ510414.1'))

    d = c.get_organism('KJ510414.1')
    print(d.nonCDS)
   # print(d)
 #   e = [] #retrive a list of organisms and iterate over them --> this will be crucial for the screener to work
''' for item in c.all_organisms().keys():
        e.append(c.get_organism(item))
    for x in e:
        for y in x.proteome:
            print(x.proteome[y][0]) # product
            print(x.proteome[y][1]) # protein location
            print(x.proteome[y][2]) # sequence'''


    
