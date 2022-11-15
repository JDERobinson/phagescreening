import re
import pandas as pd


'''Function that finds all the pam locations when given a list of pam nt seqs
and the specimens sequence. The sequence location is an exact match to the a
processed version of the fasta file given. Returns counts of pams, locations, fold
repression, and sequence of organism from fasta file'''


def get_pams(pamSeqs, fastaFile=None, sequence=None):

    pam = ['AGAAG', 'AGAAT', 'AGAAA', 'GGAAG', 'AGAAC', 'GGAAA', 'AGCAT', 'AGGAG',
    'AGGAT', 'AGCAA', 'GGAAC','GGAAT', 'AGCAG','AGGAA', 'AGGAC']

    foldrep = {'AGAAG':216.7, 'AGAAT':216.2, 'AGAAA':158.1, 'GGAAG':145.2,
    'AGAAC':120.5, 'GGAAA':110.5, 'AGCAT':84.6, 'AGGAG':82.2,
    'AGGAT':64.7, 'AGCAA':53.4, 'GGAAC':51.5, 'GGAAT':47.3, 'AGCAG':42.2,
    'AGGAA':38.5, 'AGGAC':25.5}

    if fastaFile: #if given fasta file, process and pull out sequence
        with open(fastaFile, 'r') as file: #open the kampy fasta and print out ONLY the lines
                                    #w/DNA
            text = file.readlines()
            temp = text[1:]
            seq = ""
            for line in temp:
                seq+=line
            #print(seq)
    else: seq = sequence #if given sequence, just use sequence
    PamCount = {}
    PamLoc = {}
    temp = {}
    for pam in pamSeqs:
        temp[pam] = re.finditer(pam, seq) #makes iterable dictionary with pam
                                        #sequence as key and matches (as match
                                        #objects) as the values
        locs = []
        for match in temp[pam]:
            locs.append(match.span()) #gets the location of matches within seq
            PamLoc[pam] = [foldrep[pam], locs]
            if PamCount.get(pam):
                PamCount[pam] += 1
            else:
                PamCount[pam] = 1


    return PamCount, PamLoc, seq



'''cds_nums used from EHN in 2021, edited by JDR'''

'''If given an NCBI .gb file, function returns the location, tag, and product
of all the proteins in the file'''

def cds_nums(loci_file): #loci_file is string
    cds_regions = {}
    pos_len = 0#counter for genes for cds_regions dict

    with open(loci_file, 'r', encoding = 'UTF-8') as file: #open file
        text = file.read()
        text = text.split('\n     gene')[1:]

        for entry in text:
            pos = re.search('\d*\.\.\d*', entry).group() #find position in text
            #pos in wrong format to compare to PAM seqs - current: 'start..end'
            #updated: (start, end) <- want a tuple with ints instead of just str
            pos = pos.replace('..', ',').replace('\'', '').split(',') #change
                                            #format of position to [start, end]
            positions = tuple(int(pos) for pos in pos) #make [start, end] into
                                                     #ints and convert to tuple
            locusTag = re.search('/locus_tag=".*', entry).group()[1:].replace("\"", "")
            product = re.search('/product=.*', entry).group()[1:].replace("\"", "")
            cds_regions[positions] = [locusTag, product]

    return cds_regions



    '''Function to compare pam locations to protein map in hopes of finding a
    pam location within a hypothetical protein or non-coding region

    Input: cds_nums output (dict), location of pam seqs in specimen (dict)
    Output: pandas dataframe w/ locus tag, protein product, region, and pam seq'''

    def find_pams_by_region(specimenRegions, pamLocations):
        pdtsFound = []
        pamSeqs = list(pamLocations.keys()) #makes a list of pam seqeuences -- more
                                            #accessible than if stored in dictionary
        specRegs = list(specimenRegions.keys()) #same as above
        for region in specRegs:
            for pam in pamLocations:
                seqLocation=pamLocations[pam] #made variable to make easier to ref
                for i in range(len(seqLocation)):
                    if seqLocation[i][0] >= region[0] and seqLocation[i][1] <= region[1]:
                        pdtsFound.append([specimenRegions[region][0][10:], #only tag, no label
                                        specimenRegions[region][1][8:], #only pdt no label
                                        pam, region, seqLocation])
        result_df = pd.DataFrame(data=pdtsFound,
                                columns=['LocusTag','Product','PAMseq',
                                    'SpecimenRegion','PAMLocation'])
        return result_df

    '''Finds all the pam sequences by running above function and imports into csv file'''
