#!/usr/bin/env python3
# Name: Tyler Myers (tdmyers)

import pandas as pd
import numpy as np
import re
import subprocess
import argparse
import os
# link to install hgvs module: http://hgvs.readthedocs.io/en/master/installation.html  
import hgvs.parser

#NOTE: subprocess, popen(part of subprocess) for script output capture,  
######################################################################
'''
The function of this program is to parse a BRCA exchange .tsv file for
useful delimeters and parameters. The program requires MaxEntScan's perl
script to find entropy scores for the 3 and 5 prime scores. 
'''
######################################################################



#reverse compliment of dna string
def revComp(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])

# Melissa Meredith's Codon Function (7/20/17)
def codon(BRCA, HGVS_list):    
    hgvsparser = hgvs.parser.Parser()    
    codons = {  'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C', 'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C', 'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*', 'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W', 'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R', 'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R', 'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R', 'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R','ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S', 'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S', 'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R', 'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R', 'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G', 'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G', 'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', 'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}
    codon_counter = 0
        
    # for v in range(1,len(HGVS_list)):
    #     counter+=1
    #     print("counter", counter)
    #     variant = HGVS_list[v]
    #     # run Mutalyzer here...if frame shift...
    #     var = hgvsparser.parse_hgvs_variant(variant)
    #     # the locations are integers in referance to gene start location
    #     var_Cstart_loc = var.posedit.pos.start.base
    #     var_Cend_loc = var.posedit.pos.end.base
    #     print("\tvariant", variant)
    #     print("\tHGVS:\tvar_start_loc %d \n\t var_end_loc %d" % (var_Cstart_loc, var_Cend_loc))
    #     print("var.posedit.edit", var.posedit.edit)

    for bp in range(0, len(BRCA)):    

        current_codon = ""
        protein_sequence = str() 
        codon_counter+=1
        current_codon=current_codon+BRCA[bp]
        print("current codon",( bp, current_codon))

        if (codon_counter>0 and codon_counter%3==0):
            
            if(current_codon in codons and codons[current_codon]=="*"):
                protein_sequence=protein_sequence+codons[current_codon]
                stopCodon = (bp + var_Cend_loc)
                print("\t\t(stop protein seq:", protein_sequence)
                break
            elif(current_codon in codons and codons[current_codon]!="*"):
                protein_sequence=protein_sequence+codons[current_codon]
                current_codon = ""
                codon_counter=0
            else:
                protein_sequence=protein_sequence+"!!!"
                break

    
    print("len proteinseq:\t", len(protein_sequence))
#!!!!!!!!!Return the location of the first stop codon!!!!!!!!!!!!!!!!
    if (len(protein_sequence)>0) and (protein_sequence[0]=='M' and protein_sequence[len(protein_sequence)-1]=="*"):
        print("\tcomplete coding sequence, stop codon: {}".format(protein_sequence.index('*')))
        del protein_sequence
    elif(protein_sequence[0]=='M' and protein_sequence[len(protein_sequence)-1]!="*"):
        print("\t3'-partial coding sequence")
        del protein_sequence
    elif(protein_sequence[0]!='M' and protein_sequence[len(protein_sequence)-1]=="*"):
        print("stop Codon:", stopCodon)
        print("\t 5'-partial coding sequence, no start codon (ATG: methionine) present after mutation.")
        del protein_sequence
    elif('!!!' in protein_sequence) and (protein_sequence[len(protein_sequence)-1]!="*"):
        print("\tframe-shift in coding sequence. {} amino acids after mutation location.".format(len(protein_sequence)-1))
        del protein_sequence
    else:
        print("\tinternal partial coding sequence.")
        del protein_sequence


"""This class turns all of the lines into lists, separated by tabs. the np.vstack is used
to create a matrix. this will be used later for selecting the column desired.dimesnion is
len(header)x Number of variants"""

class brcaParse:
    def __init__(self, inFile):
        f = open(inFile, "r")
        rawLabels = f.readline()
        labels = rawLabels.split("\t")

        #imports the sequence for BRCA1/2 into the program for use.
        self.BRCA1hg38Seq = open("brca1_hg38.txt", "r").read()
        self.BRCA1hg38Start = 43000000

        self.BRCA2hg38Seq = open("brca2_hg38.txt", "r").read()
        self.BRCA2hg38Start = 32300000

        self.BRCA1 = {"hg38": {"start": 43000000, "sequence": open("brca1_hg38.txt", "r").read()},
                      "hg19": {"start": 41100000, "sequence": open("brca1_hg19.txt", "r").read()}}
        self.BRCA2 = {"hg38": {"start": 32300000, "sequence": open("brca2_hg38.txt", "r").read()},
                      "hg19": {"start": 32800000, "sequence": open("brca2_hg19.txt", "r").read()}}

        #remane for better understanding!!!!
        # finds column in matrix associated with header label
        a = labels.index("Alt")
        b = labels.index("Ref")
        c = labels.index("Pos")
        d = labels.index("Clinical_significance_ENIGMA")
        e = labels.index("Gene_Symbol")
        g = labels.index("id")
        h = labels.index("Reference_Sequence")
        hgvs_cDNA = labels.index("HGVS_cDNA")

        buildMat = []
        buildMat.append(labels)

        # for loop used to append a complete data set associated with an id number. Later put into matrix/array.
        #A is all of the input data as a matrix. matOut is the selected columns as lists.
        for lines in f:
            #l = f.readline()
            if len(lines.split("\t")) == len(labels):
                buildMat.append(lines.split('\t'))
            #print(len(lines.split("\t")),len(labels))
        L = np.vstack(buildMat)
        self.A = L
        

        # use column definition to get list of a, b, c, and d. Now this is modification, reference, position and significance.
        self.Alt = brcaParse.column(self.A, a)  # positon of mutation column
        self.Ref = brcaParse.column(self.A, b)  # reference sequence column
        self.Pos = brcaParse.column(self.A, c)  # alteration of sequence column
        self.Sig = brcaParse.column(self.A, d)  # Clinical significance column
        self.Gene = brcaParse.column(self.A, e)  # BRCA1/2 column for direction of sequence
        self.id = brcaParse.column(self.A, g) #id number column
        self.ExonRef = brcaParse.column(self.A, h) #dicates what exon starts and stops should be used..reference sequence
        self.hgvs_cDNA = brcaParse.column(self.A, hgvs_cDNA) # genomic HGVS nomenclature
        self.matOut = np.vstack((self.Alt, self.Ref, self.Pos, self.Sig, self.Gene))


    #gets the desired data and puts to object
    def getDat(self):
        return (self.A, self.matOut)

    #Creates a column from the matrix of data.
    def column(matrix, i):
        return [row[i] for row in matrix]

    #Parses through the string to make all the new variant strings. If the variant is an indel
    #rather than a SNP, the definition will make more iterations to account for the new modifications.
    #tempSeq is the variant sequence, np.amax gets the maximum maxentscan score.
    def maxEntForm(self,output):
        f = open(output, 'w')
        f.write("id\tGene\tSignificance\tSpliceSite\t5'Max\t5'Ref\t3'Max\t3'Ref\tupscore\tdownscore\tinExon\n")
        for i in range(0,len(self.Gene)):
            print("i:", i)
            if self.Gene[i] == "BRCA1":
                f.write(self.id[i] +"\t" + self.Gene[i] + "\t" + self.Sig[i] + "\t")
                loc = (int(self.Pos[i]) - int(self.BRCA1hg38Start))
                #    must add a 4th variable to output the splice site later on (MMM)
                tempSeq = self.BRCA1hg38Seq[:loc-1] + self.Alt[i] + self.BRCA1hg38Seq[loc+len(self.Ref[i])-1:]

                site, upscore, downscore = self.inSpliceSite(i, tempSeq)
                f.write("{}\t".format(site))
                if (site != "N/A"):
                    upscore = 0
                    downscore = 0
                
                lenSplice = 9

                orgSeqScore, newSeqScore = self.getSeqVar(i, loc, lenSplice, tempSeq)
                if (site != "3'"):
                    f.write(str(np.amax(newSeqScore)) + "\t" + str(orgSeqScore[newSeqScore.index(np.amax(newSeqScore))]) + "\t")
                else:
                    f.write("0"+ "\t" + "0" + "\t")
                    
                lenSplice = 23
                orgSeqScore, newSeqScore = self.getSeqVar(i, loc, lenSplice, tempSeq)
                if (site != "5'"):
                    f.write(str(np.amax(newSeqScore)) + "\t" + str(orgSeqScore[newSeqScore.index(np.amax(newSeqScore))]) + "\t")
                else:
                    f.write("0"+ "\t" + "0" + "\t")
                f.write(str(upscore) + "\t" + str(downscore) + "\t")
                if(site != "N/A"):
                    f.write("1\n")
                else:
                    f.write("0\n")
               
        for i in range(0,len(self.Gene)):
            if self.Gene[i] == "BRCA2":
                f.write(self.id[i] +"\t" + self.Gene[i] + "\t" + self.Sig[i] + "\t")
                loc = (int(self.Pos[i]) - int(self.BRCA2hg38Start))
                tempSeq = self.BRCA2hg38Seq[:loc-1] + self.Alt[i] + self.BRCA2hg38Seq[loc+len(self.Ref[i])-1:]
                site, upscore, downscore = self.inSpliceSite(i, tempSeq)
                f.write("{}\t".format(site))
                if (site != "N/A"):
                    upscore = 0
                    downscore = 0
                    
                lenSplice = 9
                orgSeqScore, newSeqScore = self.getSeqVar(i, loc, lenSplice, tempSeq)
                if (site != "3'"):
                    f.write(str(np.amax(newSeqScore)) + "\t" + str(orgSeqScore[newSeqScore.index(np.amax(newSeqScore))]) + "\t")
                else:
                    f.write("0"+ "\t" + "0" + "\t")

                lenSplice = 23
                orgSeqScore, newSeqScore = self.getSeqVar(i, loc, lenSplice, tempSeq)
                if (site != "5'"):
                    f.write(str(np.amax(newSeqScore)) + "\t" + str(orgSeqScore[newSeqScore.index(np.amax(newSeqScore))]) + "\t")
                else:
                    f.write("0"+ "\t" + "0" + "\t")
                f.write(str(upscore) + "\t" + str(downscore) + "\t")
                if(site != "N/A"):
                    f.write("1\n")
                else:
                    f.write("0\n")

    def getEntScore(self,seq):
        temporary = open("temp", "w")
        temporary.write(seq)
        var = "temp"
        temporary.close()#must close file before using the subprocess
        if len(seq) == 9:
            pipe = subprocess.Popen(["perl","score5.pl", var], stdout=subprocess.PIPE)
        if len(seq) == 23:
            pipe = subprocess.Popen(["perl","score3.pl", var], stdout=subprocess.PIPE)
        result = pipe.stdout.read()
        entScore = re.findall("[+-]?\d+(?:\.\d+)?", str(result))
        return(float(entScore[0]))

    def getSeqVar(self, i, loc, lenSplice, tempSeq):
        newSeqScore = []
        orgSeqScore = []
        for j in range(0,lenSplice- 1 + (len(self.Ref[i]))):
            n = loc-lenSplice+j
            o = loc+j
            
        if self.Gene[i] =="BRCA1":
            newSeq = revComp(tempSeq[n:o])
            newSeqScore.append(self.getEntScore(newSeq))
            orgSeq = revComp(self.BRCA1hg38Seq[n:o])
            orgSeqScore.append(self.getEntScore(orgSeq))

            
        if self.Gene[i] =="BRCA2":
            newSeq = tempSeq[n:o]
            newSeqScore.append(self.getEntScore(newSeq))
            orgSeq = self.BRCA2hg38Seq[n:o]
            orgSeqScore.append(self.getEntScore(orgSeq))
        return(orgSeqScore, newSeqScore)
    
    def inSpliceSite(self, i, tempSeq):
        if self.Gene[i] == "BRCA1":
            BRCA1exons = ""
            exonStart = [43044294,43047642,43049120,43051062,43057051,43063332,43063873,43067607,43070927,
                          43074330,43076487,43082403,43090943,43094743,43095845,43097243,43099774,43104121,
                          43104867,43106455,43115725,43124016]
            exonStop = [43045802,43047703,43049194,43051117,43057135,43063373,43063951,43067695,43071238,
                        43074521,43076611,43082575,43091032,43094860,43095922,43097289,43099880,43104261,
                        43104956,43106533,43115779,43124115]
            
            print("len of self.Alt",len(self.Alt[i]))
            #all this works, butt the string is only 3682 bp's long
            if (len(self.Ref[i])!=len(self.Alt[i])) and (len(self.Alt[i])%3!=0): 
             #   print("self.Alt",self.Alt[i])
                for t in range(0,len(exonStop)):    
                    BRCA1exon = tempSeq[exonStart[t]-self.BRCA1hg38Start:exonStop[t]-self.BRCA1hg38Start]
                    BRCA1exons = BRCA1exons + BRCA1exon
                codon(revComp(BRCA1exons), self.hgvs_cDNA)

            #print("BRCA1Exon String:",len(BRCA1exons))
            
            upStream = min(exonStop, key=lambda x:abs(x-int(self.Pos[i])))
            downStream = min(exonStart, key=lambda x:abs(x-int(self.Pos[i])))
            upStreamScore, downStreamScore = self.getSpliceMaxEnt(i,upStream, downStream)
            
            for j in range(0,len(exonStart)):
                if (abs(int(self.Pos[i])-exonStart[j])<=9):
                    return("5'", upStreamScore, downStreamScore)
            for j in range(0,len(exonStop)):
                if (abs(int(self.Pos[i])-exonStop[j])<=23):
                    return("3'",upStreamScore, downStreamScore)
            else:
                return("N/A",upStreamScore, downStreamScore)
                
        if self.Gene[i] == "BRCA2":
            BRCA2exons = ""
            exonStart = [32315479,32316421,32319076,32325075,32326100,32326241,32326498,32329442,32330918,32332271,
                         32336264,32344557,32346826,32354860,32356427,32357741,32362522,32363178,32370401,32370955,
                         32376669,32379316,32379749,32380006,32394688,32396897,32398161]
            exonStop = [32315667,32316527,32319325,32325184,32326150,32326282,32326613,32329492,32331030,32333387,
                        32341196,32344653,32346896,32355288,32356609,32357929,32362693,32363533,32370557,32371100,
                        32376791,32379515,32379913,32380145,32394933,32397044,32399672]
            
            #all this works, butt the string is only 3682 bp's long
            if (len(self.Ref[i])!=len(self.Alt[i])) and (len(self.Alt[i])%3!=0): 
                print("self.Alt",self.Alt[i])
                
                for t in range(0,len(exonStop)):    
                    BRCA2exon = tempSeq[exonStart[t]-self.BRCA2hg38Start:exonStop[t]-self.BRCA2hg38Start]
                    BRCA2exons = BRCA2exons + BRCA2exon
                codon(BRCA2exons, self.hgvs_cDNA)

            #print("BRCA2Exon String:",len(BRCA2exons))
            
            upStream = min(exonStart, key=lambda x:abs(x-int(self.Pos[i])))
            downStream = min(exonStop, key=lambda x:abs(x-int(self.Pos[i])))
            upStreamScore, downStreamScore = self.getSpliceMaxEnt(i,upStream, downStream)
            
            for j in range(0,len(exonStop)):
                if (abs(int(self.Pos[i])-exonStop[j])<=9):
                    return("5'",upStreamScore, downStreamScore)
            for j in range(0,len(exonStart)):
                if (abs(int(self.Pos[i])-exonStart[j])<=23):
                    return("3'",upStreamScore, downStreamScore)
            else:
                return("N/A",upStreamScore, downStreamScore)
            
    def getSpliceMaxEnt(self, i, upStream, downStream):
        if self.Gene[i] == "BRCA1":
            loc1 = (upStream - int(self.BRCA1hg38Start))
            loc2 = (downStream - int(self.BRCA1hg38Start))
            orgScore1 = self.getEntScore(revComp(self.BRCA1hg38Seq[loc1-3:loc1+6]))
            orgScore2 = self.getEntScore(revComp(self.BRCA1hg38Seq[loc2-20:loc2+3]))
            return(orgScore1,orgScore2)

        if self.Gene[i] == "BRCA2":
            loc1 = (upStream - int(self.BRCA2hg38Start))
            loc2 = (downStream - int(self.BRCA2hg38Start))
            orgScore1 = self.getEntScore(self.BRCA2hg38Seq[loc1-3:loc1+6])
            orgScore2 = self.getEntScore(self.BRCA2hg38Seq[loc2-20:loc2+3])
            return(orgScore1,orgScore2)

                


# format for command line: pythonfile.py inputfile.tsv outputfile.tsv
def main():
    parser = argparse.ArgumentParser(description='Program to produce maxEntScan values from and input file of .tsv. output is a .tsv file as well')

    parser.add_argument('input', action="store")
    parser.add_argument('output', action="store")
    args = parser.parse_args()

    inFile = args.input
    a = brcaParse(inFile)
    A, matOut = a.getDat()
    maxEnt = args.output + ".tsv"
    a.maxEntForm(maxEnt) #get acceptor splice sites 23mers UNCOMMENT FOR 23 mer!
    os.remove("temp")

# Make loops. if BRCA1 or 2, reference BRCA1 or BRCA1, then replace string with alt(eration).

if __name__ == "__main__":
    main()
