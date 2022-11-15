
#path to file: C:\Users\Jason Dean Robinson\OneDrive\Desktop\phage lab\MycobacteriumSensor

'''Screen Kampy for PAM sequences'''

import re
from Bio import Entrez #provides tools for searching NCBI
from Bio import SeqIO #allows us to switch between sequence file types
import pandas as pd


pam = ['AGAAG', 'AGAAT', 'AGAAA', 'GGAAG', 'AGAAC', 'GGAAA', 'AGCAT', 'AGGAG',
'AGGAT', 'AGCAA', 'GGAAC','GGAAT', 'AGCAG','AGGAA', 'AGGAC']

with open('Kampy.fasta', 'r') as file: #looking at initial kampy file
    kampy = file.read()

'''Function that finds all the pam locations when given a list of pam nt seqs
and the specimens sequence. The sequence location is an exact match to the a
processed version of the fasta file given. Returns counts of pams, locations,
and sequence of organism from fasta file'''

'''210903 JDR - changed below function to take fasta file, so user does not have
to process file themselves'''

'''210906 JDR - Finished find_pams_by_region() ensuring it returns a readable
dataframe that is sorted by PAM sequence (alphabetical order). Sorting it will
help with visualization and give me the opportunity to color spreadsheet by PAM
if needed.'''

def get_pams(pamSeqs, fastaFile=None, sequence=None):

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

kamPamCount, kamPamLoc, kamSeq = get_pams(pam, "Kampy.fasta")

#print(kamPamCount)
#print(kamPamLoc)

#print(len(kampPam))

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

regions = cds_nums('Kampy_sequence.gb')
#print(regions)
#print(regions)

'''Saha gave me a region to look for pam seqeuences, so i need to pull out that
region from the text file and search it'''

#region is last protein in Kampy -- 50043:50162


prot88 = kamSeq[50043:50162] #pull out protein 88
#print(prot88)

kam88Count, kam88Loc, prot88seq = get_pams(pam, sequence = prot88)

#print('Protein 88 locations: ', kam88Loc)
#print()

'''There were zero PAM seqs in protein 88, so I am going to screen the rest of
Kampy genome for PAM seqs'''

#print('Number of PAM seqs per PAM: ', kamPamCount)
#print('Locations of PAM seqs by PAM: ', kamPamLoc)
#print('Kampy genome regions: ', regions)

'''EDIT 210902 JDR: I changed the get_pams() and cds_nums() functions so that
the sequence region output is a tuple with the start and end of of the seq as ints
ex. a gene that is from 1..4 nt is written (1, 4)

This will help when comparing Kampy regions to PAM sequence locations. The goal
is to find PAM sequences in non-coding regions or hypothetical proteins.'''

#for pam in kamPamLoc:
    #for i in range(len(kamPamLoc[pam])):
        #print(kamPamLoc[pam][i])

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

data=find_pams_by_region(regions, kamPamLoc)
#print(data)
#with open(dir, 'w') as file:
    #data.to_csv(dir)

'''There are 150 PAM positions within hypothetical proteins. TO DO: Look at df
to figure out which PAM seq within a protein we should knock out -- consult Saha
on this matter'''

'''210916 JDR'''

'''210924 JDR - just moved some code around so that the original data set is
being edited (add FoldRepression, etc) instead of the just the hypothetical PAMS
db. Also created a non-hypo GFP-length db.'''

'''Add FoldRepression to data sheet'''
def fold(pam):
    foldrep = {'AGAAG':216.7, 'AGAAT':216.2, 'AGAAA':158.1, 'GGAAG':145.2,
    'AGAAC':120.5, 'GGAAA':110.5, 'AGCAT':84.6, 'AGGAG':82.2,
    'AGGAT':64.7, 'AGCAA':53.4, 'GGAAC':51.5, 'GGAAT':47.3, 'AGCAG':42.2,
    'AGGAA':38.5, 'AGGAC':25.5}
    repression = foldrep[pam]
    return repression

result = data['PAMseq'].apply(fold)
#print(result)
data['FoldRepression'] = result
#print(data.dtypes)

'''Split SpecimenRegion into RegionStart and RegionEnd'''
a = '(123, 456)'
a = a.split(',')
first = int(a[0].replace(')', "").replace('(', "").replace(" ", ""))
second = int(a[1].replace(')', "").replace('(', "").replace(" ", ""))
#print(first, second)

pamDat = data['SpecimenRegion'].astype('str')
pamDat = pamDat.str.split(',')
#print(pamDat)
start = pamDat.apply(lambda x: int(x[0].replace(')', "").replace('(', "").replace(" ", "")))
#print(start)
end = pamDat.apply(lambda x: int(x[1].replace(')', "").replace('(', "").replace(" ", "")))
#print(end)

data['RegionStart'] = start
data['RegionEnd'] = end
#print(data.head)
'''Drop any data with FoldRepression less than Hatfull threshold: 25'''

data=data.sort_values('PAMseq')
data = data[data['FoldRepression'] >= 25 ]
dir='C:\\Users\\Jason Dean Robinson\\OneDrive\\Desktop\\phage lab\\MycobacteriumSensor\\PAMSeqs.csv'
#data.to_csv(dir)

'''Processes dataframe to pull out only the lines that are hypothetical
proteins'''

hypoPams=data[(data['Product'] == 'hypothetical protein')]
#print(hypoPams)
#print(type(hypoPams))
#print(hypoPams.shape)
hypoPams=hypoPams.sort_values('PAMseq')
dir='C:\\Users\\Jason Dean Robinson\\OneDrive\\Desktop\\phage lab\\MycobacteriumSensor\\hypoPams.csv'
#hypoPams.to_csv(dir)

newdir='C:\\Users\\Jason Dean Robinson\\OneDrive\\Desktop\\phage lab\\MycobacteriumSensor\\nohypoPams.csv'
noHypo= data[(data['Product']!='hypothetical protein')]
#noHypo.to_csv(newdir)

'''Now drop any rows whose bp sequence is less than 760 bp or greater than 1000
and create spreadsheet with them.'''

#hypothetical regions
hypoPams['TotalLen'] = hypoPams['RegionEnd'] - hypoPams['RegionStart']
#print(hypoPams['TotalLen'])

gfpPams = hypoPams[(hypoPams['TotalLen'] >= 750) & (hypoPams['TotalLen'] <= 1000)]
gfpPams = gfpPams.sort_values('FoldRepression', ascending=False)

#print(gfpPams)
dir = 'C:\\Users\\Jason Dean Robinson\\OneDrive\\Desktop\\phage lab\\MycobacteriumSensor\\HypoGfpPams.csv'
#gfpPams.to_csv(dir)

#non-hypothetical regions
noHypo['TotalLen'] = noHypo['RegionEnd'] - noHypo['RegionStart']
#print(noHypo['TotalLen'])

gfpPamsNon = noHypo[(noHypo['TotalLen'] >= 750) & (noHypo['TotalLen'] <= 1000)]
gfpPamsNon = gfpPamsNon.sort_values('FoldRepression', ascending=False)

#print(gfpPamsNon)
dir = 'C:\\Users\\Jason Dean Robinson\\OneDrive\\Desktop\\phage lab\\MycobacteriumSensor\\gfpPams.csv'
gfpPamsNon.to_csv(dir)

'''Pull out the regions we will be replacing into 2 variables. Later I need
to add the 100 nt 5' and 3' of gene'''

immuneRep = noHypo[(noHypo['Product'] == 'immunity repressor')]
column = immuneRep['FoldRepression']
max = column.max()
immuneRep = immuneRep[column == max] #pull out only the immunity repressor rows that have the max PAM fold repression

locIM = []
locIM = [int(immuneRep['RegionStart']), int(immuneRep['RegionEnd'])]
#print('Immune Repressor Location:', locIM)
#print(immuneRep.dtypes)
pamIM=immuneRep.iloc[0]['PAMseq'] #all the immune suppresor PAMS in IR with a max fold repression
print(pamIM, 'are the PAMS')
#there is only one PAM within the IR sequence: AGGAG

pamIM, foldrepIM, pamsLoc=immuneRep.iloc[0]['PAMseq'], immuneRep.iloc[0]['FoldRepression'], immuneRep.iloc[0]['PAMLocation']



for item in pamsLoc: #this is just getting exta information about the PAM, including location and fold repression
    if (item[0] >= locIM[0]) and item[1] <= locIM[1]:
        pamsLoc = item
print(pamIM, foldrepIM, pamsLoc)
#okay, so the df I am referencing has the PAM in one column, and it's list of locations.
#currently pamsLoc includes the pam location for entire genome, not just immunity
#repressor location, so I have to get the location of my PAM sequence WITHIN
#the immunity repressor


lysinB = gfpPamsNon[(gfpPamsNon['Product'] == 'lysin B')]
lysinB = lysinB.iloc[[0]]
locLys = [int(lysinB['RegionStart']), int(lysinB['RegionEnd'])]
#print('Lysin B Location:', locLys)

##EDIT 220329 - cleaning up some code, adding comments, fixing slicing of DNA
'''Pulling out gene from fasta with 100 nt 5' and 3 '''

with open('Kampy.fasta', 'r', encoding = 'UTF-8') as file: #open file
    next(file)
    text = file.read().replace('\n', '')
    #print(text)
    lys = text[(locLys[0]-99):(locLys[1]+101)]
    #print(lys)
    immune = text[(locIM[0]-99):(locIM[1]+101)]
    #print()
    print(immune)
    #print('LysinB :', len(lys))
    print('Immunity Repressor Len:', len(immune))

'''220120: Now that I know the location of the PAM sequence within IR, can pull
out the sgRNA sequence'''

pamIRloc = re.search('AGGAG', immune).span()
print(pamIRloc)
sgRNA = immune[pamIRloc[1]:pamIRloc[1]+21].replace('\n', '')
print(sgRNA)
print(len(sgRNA))

#NOTE: Need to check if the sequences are correct, perhaps check to see if the
#lengths are what is expected by comparing to TotalLength column for that gene
print()

#get just the homolgous sections up and downstream of ea gene and put into a dict
#downstream is 0:100 from 5' to 3'
#upstream is 100 from the end, to the end 5' to 3'

upLys = lys[:101].replace('\n', '')
downLys = lys[len(lys)-101:].replace('\n', '')

upIm = immune[:101].replace('\n', '')
downIm = immune[len(immune)-101:].replace('\n', '')

lysHom = {'Gene': 'LysinB', 'upstream':upLys, 'downstream':downLys}
immuneHom = {'Gene':'ImmunityRepressor', 'upstream':upIm, 'downstream':downIm}
print('Lysin B', lysHom)
print()
print('Immune Repressor', immuneHom, end = '\n')
print()

upstream=['UPSTREAM\nhomologous oligo 1: '+downIm[50:]+
'\nhomologous oligo 2: '+downIm[0:50]+
'\noverlap between oligo 1 and 2: '+downIm[50:75]]

downstream = ['DOWNSTREAM\n'+
'homologous oligo 1: '+upIm[0:50]+
'\nhomologous oligo 2: '+upIm[50:]+
'\noverlap between oligo 1 and 2: '+upIm[25:50]]

'''211201 I'm just going to add everything to a text file so I don't have to keep
running this file to get the oligonucleotides'''

with open('Oligo.txt', 'w+') as file:
    file.write('Location: '+str(locIM)+'\n\n')
    file.write('ImmunityRepressor Gene:\n'+immune)
    file.write('\n'*2)
    file.writelines(upstream)
    file.write('\n'*2)
    file.writelines(downstream)
    file.write('\n'*2)
    file.writelines(['PAM sequence:\n', pamIM,'\n'*2, 'Fold Repression:\n', str(foldrepIM)])
    file.write('\n'*2)
    file.writelines(['PAM location:\n',str(pamIRloc),'\n'*2,'sgRNA sequence:\n',sgRNA])

    #file.write()
