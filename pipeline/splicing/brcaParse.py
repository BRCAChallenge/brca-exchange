#!/usr/bin/env python3
# Name: Tyler Myers (tdmyers)

import pandas as pd
import numpy as np
import re
import subprocess
import argparse
import os
import exonDict

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

#Calls the perl scripts to get the MaxEntScore of the dna sequence
def getEntScore(seq):
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
        g = labels.index("HGVS_cDNA")
        h = labels.index("Reference_Sequence")

        buildMat = []
        buildMat.append(labels)

        # for loop used to append a complete data set associated with an id number. Later put into matrix/array.
        #A is all of the input data as a matrix. matOut is the selected columns as lists.
        for lines in f:
            #l = f.readline()
            if len(lines.split("\t")) == len(labels):
                buildMat.append(lines.split('\t'))
        L = np.vstack(buildMat)
        self.A = L
        

        # use column definition to get list of a, b, c, and d. Now this is modification, reference, position and significance.
        self.Alt = brcaParse.column(self.A, a)  # positon of mutation column
        self.Ref = brcaParse.column(self.A, b)  # reference sequence column
        self.Pos = brcaParse.column(self.A, c)  # alteration of sequence column
        self.Sig = brcaParse.column(self.A, d)  # Clinical significance column
        self.Gene = brcaParse.column(self.A, e)  # BRCA1/2 column for direction of sequence
        self.id = brcaParse.column(self.A, g) #id HGVS_cDNA column
        self.ExonRef = brcaParse.column(self.A, h) #dicates what exon starts and stops should be used..reference sequence
        
        self.matOut = np.vstack((self.Alt, self.Ref, self.Pos, self.Sig, self.Gene))

    #Creates a column from the matrix of data.
    def column(matrix, i):
        return [row[i] for row in matrix]
    

    #Parses through the string to make all the new variant strings. If the variant is an indel
    #rather than a SNP, the definition will make more iterations to account for the new modifications.
    #tempSeq is the variant sequence, np.amax gets the maximum maxentscan score.
    def maxEntForm(self,output):
        f = open(output, 'w')
        f.write("HGVS_cDNA\tGene\tSignificance\tSpliceSite\t5'Alt\t5'AltZScore\t5'Ref\t5'RefZScore\t3'Alt\t3'AltZScore\t3'Ref\t3'RefZScore\tupscore\tdownscore"+
                "\tMax5'deNovo\tMax5'deNovoZScore\t5'Ref\t5'ZScore\tmaxScoreLoc\tinSpliceSite\texonLoc\tpathProb\n")
        for i in range(0,len(self.Gene)):
            if self.Gene[i] == "BRCA1":
                f.write(self.id[i] +"\t" + self.Gene[i] + "\t" + self.Sig[i] + "\t")
                loc = (int(self.Pos[i]) - int(self.BRCA1hg38Start))
                site, upscore, downscore, exonLoc = self.inSpliceSite(i)
                exLoc = int(exonLoc) - int(self.BRCA1hg38Start)#normalized location for seuqence
                f.write("{}\t".format(site))

                if (site != "N/A" or site !="inExon"):
                    upscore = 0
                    downscore = 0
                    
                lenSplice = 9
                tempSeq = self.BRCA1hg38Seq[:loc-1] + self.Alt[i] + self.BRCA1hg38Seq[loc+len(self.Ref[i])-1:]
                if (site == "5'"):
                    altScore = getEntScore(revComp(tempSeq[exLoc-3:exLoc+6]))#score for the splice site with the mutation
                    oriScore = getEntScore(revComp(self.BRCA1hg38Seq[exLoc-3:exLoc+6]))#score for unaltered cDNA
                    

                    f.write(str(altScore) + "\t" + str(self.getZScore(altScore,site))+
                            "\t"+ str(oriScore) + "\t"+ str(self.getZScore(oriScore,site))+"\t")
                    
                    pathProb = self.pathProb(self.getZScore(oriScore,site),self.getZScore(altScore,site),site,i)
                else:
                    f.write("0"+ "\t" + "0" + "\t"+"0"+ "\t" + "0" + "\t")

                lenSplice = 23

                if (site == "3'"):
                    altScore = getEntScore(revComp(tempSeq[exLoc-20:exLoc+3]))#score for the splice site with the mutation
                    oriScore = getEntScore(revComp(self.BRCA1hg38Seq[exLoc-20:exLoc+3]))#score for unaltered cDNA
                    
                    f.write(str(altScore) + "\t" + str(self.getZScore(altScore,site))+
                            "\t"+ str(oriScore) + "\t"+ str(self.getZScore(oriScore,site))+"\t")
                    
                    pathProb = self.pathProb(self.getZScore(oriScore,site),self.getZScore(altScore,site),site,i)
                else:
                    f.write("0"+ "\t" + "0" + "\t"+"0"+ "\t" + "0" + "\t")
                f.write(str(upscore) + "\t" + str(downscore) + "\t")

                lenSplice = 9
                if(site == "N/A" or site =="inExon"):
                    orgSeqScore, newSeqScore = self.getSeqVar(i, loc, lenSplice, tempSeq)
                    pathProb = self.pathProb(self.getZScore(orgSeqScore[newSeqScore.index(np.amax(newSeqScore))],site),self.getZScore(np.amax(newSeqScore),site),site,i)

                    f.write(str(np.amax(newSeqScore)) + "\t" + str(self.getZScore(np.amax(newSeqScore),site))+
                            "\t"+ str(orgSeqScore[newSeqScore.index(np.amax(newSeqScore))]) + "\t"+ str(self.getZScore(orgSeqScore[newSeqScore.index(np.amax(newSeqScore))],site))+"\t")
                    index = int(self.Pos[i]) + newSeqScore.index(np.amax(newSeqScore))
                    f.write(str(index)+"\t")


                if(site != "N/A" and site !="inExon"):
                    f.write("0\t0\t0\t0\t0\t1\t" +str(exonLoc)+"\t")
                else:
                    f.write("0\t"+str(exonLoc)+"\t")
                f.write(str(pathProb) +"\n")
                
               
        for i in range(0,len(self.Gene)):
            if self.Gene[i] == "BRCA2":
                f.write(self.id[i] +"\t" + self.Gene[i] + "\t" + self.Sig[i] + "\t")
                loc = (int(self.Pos[i]) - int(self.BRCA2hg38Start))
                site, upscore, downscore, exonLoc = self.inSpliceSite(i)
                exLoc = int(exonLoc) - int(self.BRCA2hg38Start)#normalized location for seuqence
                f.write("{}\t".format(site))

                if (site != "N/A" or site !="inExon"):
                    upscore = 0
                    downscore = 0
                    
                lenSplice = 9
                tempSeq = self.BRCA1hg38Seq[:loc-1] + self.Alt[i] + self.BRCA1hg38Seq[loc+len(self.Ref[i])-1:]
                if (site == "5'"):
                    altScore = getEntScore(tempSeq[exLoc-3:exLoc+6])#score for the splice site with the mutation
                    oriScore = getEntScore(self.BRCA2hg38Seq[exLoc-3:exLoc+6])#score for unaltered cDNA
                    

                    f.write(str(altScore) + "\t" + str(self.getZScore(altScore,site))+
                            "\t"+ str(oriScore) + "\t"+ str(self.getZScore(oriScore,site))+"\t")
                    
                    pathProb = self.pathProb(self.getZScore(oriScore,site),self.getZScore(altScore,site),site,i)
                else:
                    f.write("0"+ "\t" + "0" + "\t"+"0"+ "\t" + "0" + "\t")

                lenSplice = 23
                pathProb = self.pathProb(self.getZScore(orgSeqScore[newSeqScore.index(np.amax(newSeqScore))],site),self.getZScore(np.amax(newSeqScore),site),site,i)

                if (site == "3'"):
                    altScore = getEntScore(tempSeq[exLoc-20:exLoc+3])#score for the splice site with the mutation
                    oriScore = getEntScore(self.BRCA2hg38Seq[exLoc-20:exLoc+3])#score for unaltered cDNA
                    
                    f.write(str(altScore) + "\t" + str(self.getZScore(altScore,site))+
                            "\t"+ str(oriScore) + "\t"+ str(self.getZScore(oriScore,site))+"\t")
                else:
                    f.write("0"+ "\t" + "0" + "\t"+"0"+ "\t" + "0" + "\t")
                f.write(str(upscore) + "\t" + str(downscore) + "\t")

                lenSplice = 9
                if(site == "N/A" or site =="inExon"):
                    orgSeqScore, newSeqScore = self.getSeqVar(i, loc, lenSplice, tempSeq)
                    pathProb = self.pathProb(self.getZScore(orgSeqScore[newSeqScore.index(np.amax(newSeqScore))],site),self.getZScore(np.amax(newSeqScore),site),site,i)

                    f.write(str(np.amax(newSeqScore)) + "\t" + str(self.getZScore(np.amax(newSeqScore),site))+
                            "\t"+ str(orgSeqScore[newSeqScore.index(np.amax(newSeqScore))]) + "\t"+ str(self.getZScore(orgSeqScore[newSeqScore.index(np.amax(newSeqScore))],site))+"\t")
                    index = int(self.Pos[i]) + newSeqScore.index(np.amax(newSeqScore))
                    f.write(str(index)+"\t")

                if(site != "N/A" and site !="inExon"):
                    f.write("0\t0\t0\t0\t0\t1\t" +str(exonLoc)+"\t")
                else:
                    f.write("0\t"+str(exonLoc)+"\t")
                f.write(str(pathProb) +"\n")
                
  

    def getSeqVar(self, i, loc, lenSplice, tempSeq):
        newSeqScore = []
        orgSeqScore = []
        for j in range(0,lenSplice- 1 + (len(self.Ref[i]))):
            n = loc-lenSplice+j
            o = loc+j
            if self.Gene[i] =="BRCA1":
                newSeq = revComp(tempSeq[n:o])
                newSeqScore.append(getEntScore(newSeq))
                orgSeq = revComp(self.BRCA1hg38Seq[n:o])
                orgSeqScore.append(getEntScore(orgSeq))

            
            if self.Gene[i] =="BRCA2":
                newSeq = tempSeq[n:o]
                newSeqScore.append(getEntScore(newSeq))
                orgSeq = self.BRCA2hg38Seq[n:o]
                orgSeqScore.append(getEntScore(orgSeq))
        return(orgSeqScore, newSeqScore)
    
    def inSpliceSite(self, i):
        if self.Gene[i] == "BRCA1":
            exonStart = exonDict.exonStarts.get(self.ExonRef[i])
            exonStop = exonDict.exonStops.get(self.ExonRef[i])
            upStream = min(exonStop, key=lambda x:abs(x-int(self.Pos[i])))
            exonStop.remove(upStream)
            downStream = min(exonStop, key=lambda x:abs(x-int(self.Pos[i])))
            exonStop.append(upStream)
            upStreamScore, downStreamScore = self.getSpliceMaxEnt(i,upStream, downStream)

            
            for j in range(0,len(exonStart)):
                if (abs(int(self.Pos[i])-exonStart[j])<=9):
                    exonLoc = exonStart[j]
                    return("5'", upStreamScore, downStreamScore,exonLoc)

            for j in range(0,len(exonStop)):
                if (abs(int(self.Pos[i])-exonStop[j])<=23):
                    exonLoc = exonStop[j]
                    return("3'",upStreamScore, downStreamScore,exonLoc)
                
            for j in range(0,len(exonStart)):
                if (int(self.Pos[i])>=exonStart[j] and int(self.Pos[i])<=exonStop[j]):
                    return("inExon",upStreamScore, downStreamScore,0)
            else:
                return("N/A",upStreamScore, downStreamScore,0)
                
        if self.Gene[i] == "BRCA2":
            exonStart = exonDict.exonStarts.get(self.ExonRef[i])
            exonStop = exonDict.exonStops.get(self.ExonRef[i])
            
            upStream = min(exonStart, key=lambda x:abs(x-int(self.Pos[i])))
            exonStart.remove(upStream)
            downStream = min(exonStart, key=lambda x:abs(x-int(self.Pos[i])))
            exonStart.append(upStream)
            upStreamScore, downStreamScore = self.getSpliceMaxEnt(i,upStream, downStream)
           
            for j in range(0,len(exonStop)):
                if (abs(int(self.Pos[i])-exonStop[j])<=9):
                    exonLoc = exonStart[j]
                    return("5'",upStreamScore, downStreamScore, exonLoc)
                
            for j in range(0,len(exonStart)):
                if (abs(int(self.Pos[i])-exonStart[j])<=23):
                    exonLoc = exonStart[j]                   
                    return("3'",upStreamScore, downStreamScore, exonLoc)

            for j in range(0,len(exonStart)):
                if (int(self.Pos[i])>=exonStart[j] and int(self.Pos[i])<=exonStop[j]):
                    return("inExon",upStreamScore, downStreamScore,0)
            else:
                return("N/A",upStreamScore, downStreamScore,0)
            
    def getSpliceMaxEnt(self, i, upStream, downStream):
        if self.Gene[i] == "BRCA1":

            loc1 = (upStream - int(self.BRCA1hg38Start))
            loc2 = (downStream - int(self.BRCA1hg38Start))
            orgScore1 = getEntScore(revComp(self.BRCA1hg38Seq[loc1-3:loc1+6]))
            orgScore2 = getEntScore(revComp(self.BRCA1hg38Seq[loc2-3:loc2+6]))
            return(orgScore1,orgScore2)

        if self.Gene[i] == "BRCA2":
            
            loc1 = (upStream - int(self.BRCA2hg38Start))
            loc2 = (downStream - int(self.BRCA2hg38Start))
            orgScore1 = getEntScore(self.BRCA2hg38Seq[loc1-3:loc1+6])
            orgScore2 = getEntScore(self.BRCA2hg38Seq[loc2-3:loc2+6])
            return(orgScore1,orgScore2)

    def getZScore(self,entScore,site):
        #set up parse txt file
        
        donorstd = 2.16#240#5334713458
        donormean = 8.22#2291666666667

        acceptorstd = 2.61#74666704340925
        acceptormean =8.13#87499999999999
        if(site =="5'"):
            score = (entScore-donormean)/donorstd
             
        if (site =="3'"):
            score = (entScore-acceptormean)/acceptorstd

        if(site == "N/A" or site =="inExon"):
            score = (entScore-donormean)/donorstd
        return(score)

    def pathProb(self,zScoreRef, zScoreAlt, site,i):
        pathProb = 0
        if (site == "3'"):
            if (zScoreAlt > zScoreRef):
                pass
            if (zScoreRef < -1.5 and zScoreAlt-zScoreRef>0.5):
                pathProb = 0.97
            if (zScoreAlt > 0.5):
                pathProb = 0
            if (zScoreAlt >= -1.5 and zScoreAlt <= 0.5):
                pathProb = 0.34
            if (zScoreAlt < -1.5):
                pathProb = 0.97
        if (site == "5'"):
            if (zScoreAlt > zScoreRef):
                pass
            if (zScoreRef < -1.0 and zScoreAlt-zScoreRef>0.5):
                pathProb = 0.97
            if (zScoreAlt > 0):
                pathProb = 0
            if (zScoreAlt >= -2 and zScoreAlt <= 0):
                pathProb = 0.34
            if (zScoreAlt < -2):
                pathProb = 0.97
        if (site == "N/A" or site =="inExon"):
            exonStart = exonDict.exonStarts.get(self.ExonRef[i])
            exonStop = exonDict.exonStops.get(self.ExonRef[i])

            for j in range(0,len(exonStart)):
                if(int(self.Pos[i])>=exonStart[j] and int(self.Pos[i])<=exonStop[j]):
                    if zScoreAlt <-2:
                        pass
                    if zScoreAlt>= -2 and zScoreAlt<0.3:
                        pathProb = 0.3
                    if zScoreAlt>0:
                        pathProb = 0.64
#if zref(subsequent splice donor)< zScoreAlt, promote path prob by one step
        return(pathProb)
            
            


# format for command line: pythonfile.py inputfile.tsv outputfile.tsv
def main():
    parser = argparse.ArgumentParser(description='Program to produce maxEntScan values from and input file of .tsv. output is a .tsv file as well')

    parser.add_argument('input', action="store")
    parser.add_argument('output', action="store")
    args = parser.parse_args()

    inFile = args.input
    a = brcaParse(inFile)
    maxEnt = args.output + ".tsv"
    a.maxEntForm(maxEnt) #get acceptor splice sites 23mers UNCOMMENT FOR 23 mer!
    os.remove("temp")

# Make loops. if BRCA1 or 2, reference BRCA1 or BRCA1, then replace string with alt(eration).

if __name__ == "__main__":
    main()
