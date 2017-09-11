
# Directory for splice site analysis code.
brcaParse.py gives MaxEntScan scores and other positional information about mutations in BRCA1/2 genes. The program uses the .tsv output from the BRCA Exchange webpage and outputs a output.tsv with splicing information about the genetic variants. 

The command line format is:

python3 input.tsv output.tsv

#__init__:
This definition creates strings for the BRCA1 and BRCA2 sequences. self.BRCA1hg38 and self.BRCA2hg38 are the sequences used for the reference gene sequences. a, b, c, d, e, g, and h are used as the labels to strip from the .tsv file for this program. The class attributes Alt, Ref, Pos, Sig, Gene, id, and ExonRef are used throughout this program.

Alt: the alteration/mutation in the genomic sequence.
Ref: the reference character in the genomic sequence.
Pos: position of mutation or indel.
Sig: clinical significance of the mutation. 
Gene: dictates the gene(either BRCA1 or 2).
id: the id number from the BRCAExchange.
ExonRef: the cDNA reference number to find the corresponding Exons. 

#def column:
takes in a column from the .tsv and turns the column into a list

#def maxEntForm:

the file, f, is opened and assigned the name given from the command line. The name comes from the output.tsv field of the command line. the program writes the header for the fields of the output tsv. 

The for loop over the length of self.Gene keeps track of the individual variant id. The position and id are written to tsv. 

in the for loop, the program checks for BRCA1 or BRCA2 genes to differentiate between forward and reverse strand. The program then checks if the mutation is in a splice site as well as the scores of the neighboring splice sites. We are then informed by the variable site, if the site is de novo, 5' or 3'. 

The lenSplice is set to 9 for the 5' splice sites and then the tempSeq is created using the alteration column of the tsv file. The self.GetSeqVar function then gets the new and old maxentscores for the mutation. 

The score is then written to the .tsv, along with the z scores and pathProbability.

The same process is then repeated with the 23 mer 3' site. 

After the 5' and 3' procedure is repeated with BRCA2 labeled genes. 

#def getSeqVar:
getSeqVar retrieves all the MaEntScores for the positions affected by the mutation.

BRCA1 is on the reverse strand so revComp of the string is created for computation.

#def inSpliceSite:
exonStart and exonStop are generated from a dictionary of exon starts and stops from the file exonDict.py. the cDNA reference number is used to select the correct starts and stops. upStreamScore and downStreamScore are the scores of the upstream and downstream splice sites. 

The following for loops detect if the mutations are within the bounds of the splice site. 

#def getSpliceMaxEnt:
getSpliceMaxEnt retrieves the maxEntScore of the upstream and down stream splice sites.

#getZScore:
given the mean, standard deviation, entScore, and site type, getZScore calculates the splice sites z score. 

#def pathProb:
following the research on BRCA1/2 pathogenicity Probability, pathProb decides the pathProb based on the z scores of the maxEntScores based on the comparison of the reference sequences and the mutated sequence. 

for "N/A" splice sites, the program verifies that the mutation is within an exon and calculates the donor splice score. 



