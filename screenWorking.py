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
                self.proteome = self.__chunk_fasta(self.genome, orgProteome)
            else: self.proteome = None
        
        def __chunk_fasta(self, cfGenome, cfProteome):
            chunks = {}

            for protein in cfProteome:
               # print(cfProteome[protein])
                cds = []
                start = cfProteome[protein][1][0]
                end = cfProteome[protein][1][1] 
                pdt = cfProteome[protein][0]
              #  print(pdt)
                cds = [pdt, (start, end), cfGenome[start:end]]
                chunks[protein] = cds 
            
            return chunks
        
        def __str__(self): ##summary statement of an organism, prints ID, length of genome, and number of genes + tags
            if self.proteome:
              #  s = (f'Reference ID: {self.id}\nGenome length: {len(self.genome)}\nNumber of proteins: {len(self.proteome)}',
                    # f'\nProteins: {[item for item in self.proteome[item][0]]}') 
           
               # return 
            
   
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

   # d = c.get_organism('KJ510414.1')
   # print(type(d.proteome))

    e = c.all_organisms()
    print(e['KJ510414.1'].proteome)