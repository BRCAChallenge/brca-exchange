#!/usr/bin/env python
"""calcMaxEntScanMeanStd: calculate the mean and stddev of the donor and acceptor maxEntScan scores
   over the indicated transcripts.  These scores can be used to convert raw maxEntscan scores into z-scores.
   Generates an output file in json format.

   Usage: calcMaxEntscanMeanStd -t NM_007294.3 NM_000059.3  -o brca.zscore.json

"""
import argparse
from Bio.Seq import Seq
import json
from MySQLdb.constants import FIELD_TYPE
import _mysql
import numpy
import os
import re
import requests
import subprocess
import tempfile

# basically a faster version of score3/score5
from maxentpy import maxent
from maxentpy import maxent_fast  # fast version could further speed up
from maxentpy.maxent_fast import load_matrix  # support preloading matrix
matrix3 = load_matrix(3)
matrix5 = load_matrix(5)

db = None


# static gene coordinate data for failed calls to genome-mysql.cse.ucsc.edu
TranscriptDataBRCA1 = {'bin': '114',
                       'exonEnds': '43045802,43047703,43049194,43051117,43057135,43063373,43063951,43067695,43071238,43074521,43076614,43082575,43091032,43094860,43095922,43097289,43099880,43104261,43104956,43106533,43115779,43124115,43125483,',
                       'exonFrames': '1,0,1,0,0,1,1,0,1,2,1,0,1,1,2,1,0,1,2,2,2,0,-1,',
                       'name': 'NM_007294.3',
                       'txStart': 43044294,
                       'exonCount': 23,
                       'cdsEndStat': 'cmpl',
                       'cdsEnd': 43124096,
                       'score': 0,
                       'name2': 'BRCA1',
                       'strand': '-',
                       'cdsStart': 43045677,
                       'cdsStartStat': 'cmpl',
                       'chrom': 'chr17',
                       'txEnd': 43125483,
                       'exonStarts': '43044294,43047642,43049120,43051062,43057051,43063332,43063873,43067607,43070927,43074330,43076487,43082403,43090943,43091434,43095845,43097243,43099774,43104121,43104867,43106455,43115725,43124016,43125270,'}

TranscriptDataBRCA2 = {'bin': '103',
                       'exonEnds': '32315667,32316527,32319325,32325184,32326150,32326282,32326613,32329492,32331030,32333387,32341196,32344653,32346896,32355288,32356609,32357929,32362693,32363533,32370557,32371100,32376791,32379515,32379913,32380145,32394933,32397044,32399672,',
                       'exonFrames': '-1,0,1,1,2,1,0,1,0,1,1,1,1,2,1,0,2,2,0,0,1,0,1,0,1,0,0,',
                       'name': 'NM_000059.3',
                       'txStart': 32315479,
                       'exonCount': 27,
                       'cdsEndStat': 'cmpl',
                       'cdsEnd': 32398770,
                       'score': 0,
                       'name2': 'BRCA2',
                       'strand': '+',
                       'cdsStart': 32316460,
                       'cdsStartStat': 'cmpl',
                       'chrom': 'chr13',
                       'txEnd': 32399672,
                       'exonStarts': '32315479,32316421,32319076,32325075,32326100,32326241,32326498,32329442,32330918,32332271,32336264,32344557,32346826,32354860,32356427,32357741,32362522,32363178,32370401,32370955,32376669,32379316,32379749,32380006,32394688,32396897,32398161,'}


def fetch_gene_coordinates(transcript_name):
    """Query the indicated transcripts from the genome browser hg38 refseq table"""
    try:
        global db # db is global to prevent reconnecting.
        if db is None:
            print 'connect'
        conv= { FIELD_TYPE.LONG: int }
        db = _mysql.connect(host='genome-mysql.cse.ucsc.edu',user='genome',passwd='',db="hg38",conv=conv)
        db.query("""SELECT * FROM ncbiRefSeq WHERE name = '%s'""" % transcript_name)
        r = db.use_result().fetch_row(how=1,maxrows=0)
        if len(r)>1:
            pass
        else:
            return r[0]
    except IndexError as e:
        if transcript_name == 'NM_007294.3':
            return TranscriptDataBRCA1
        elif transcript_name == 'NM_000059.3':
            return TranscriptDataBRCA2



def runMaxEntScan(sequence, donor=False, usePerl=False):
    """Run maxEntScan on the indicated sequence.  Run score5.pl on candidate donor sequences, 
       score3.pl on candidate acceptor sequences.  Return the score"""
    if usePerl:
        (fd, tmpfile) = tempfile.mkstemp()
        fp = os.fdopen(fd, "w")
        fp.write(sequence)
        fp.close()
        if donor:
            pipe = subprocess.Popen(["perl", os.path.join(os.path.dirname(__file__), 'score5.pl'), tmpfile], stdout=subprocess.PIPE)
        else:
            pipe = subprocess.Popen(["perl", os.path.join(os.path.dirname(__file__), 'score3.pl'), tmpfile], stdout=subprocess.PIPE)
        result = pipe.stdout.read()
        entScore = re.findall("[+-]?\d+(?:\.\d+)?", str(result))
        os.remove(tmpfile)
        return(float(entScore[0]))
    else:
        # use the fastest available maxentpy implementation
        # we round to two decimal places to match the perl script's output
        if donor:
            return round(maxent_fast.score5(str(sequence), matrix=matrix5), 2)
        else:
            return round(maxent_fast.score3(str(sequence), matrix=matrix3), 2)


def scoreSeq(chrom, strand, coordinate, donor=False, verbose=False):
    """Given the coordinate of a putative splice site (either a donor or acceptor, as indicated by
       the donor argument), return the maxEntScan score describing the splice site fitness of the genome at
       that coordinate"""
    if donor:
        if strand == '+':
            rangeStartCoord = int(coordinate) - 3 + 1
            rangeEndCoord = int(coordinate) + 6
        else:
            rangeStartCoord = int(coordinate) - 6 + 1
            rangeEndCoord = int(coordinate) + 3
    else:
        if strand == '+':
            rangeStartCoord = int(coordinate) - 20 + 1
            rangeEndCoord = int(coordinate) + 3
        else:
            rangeStartCoord = int(coordinate) - 3 + 1
            rangeEndCoord = int(coordinate) + 20
    url = "http://togows.org/api/ucsc/hg38/%s:%d-%d.fasta" % (chrom, rangeStartCoord, rangeEndCoord)
    rr = requests.get(url)
    if verbose:
        print rr.content
    lines = rr.content.split('\n')
    headerLine = lines[0]
    sequence = ""
    for ii in range(1,len(lines)):
        sequence += lines[ii]
    if strand == '-':
        sequence = str(Seq(sequence).reverse_complement())
    if verbose:
        print "Sequence to be scored:", sequence
    score = runMaxEntScan(sequence, donor=donor)
    if verbose:
        print "score:", score
    return score

def addDataForThisTranscript(donors, acceptors, transcript, verbose=False):
    """Given a transcript, score the interior splice sites of that transcript, and add the scores
       to the donor and acceptor arrays respectively."""
    if transcript['strand'] == '+':
        exonStarts = transcript['exonStarts'].split(',')
        exonEnds = transcript['exonEnds'].split(',')
    else:
        exonStarts = list(reversed(transcript['exonEnds'].split(',')))
        exonEnds = list(reversed(transcript['exonStarts'].split(',')))
    for ii in range(transcript['exonCount']):
        if ii == 0:
            donors = numpy.append(donors, [scoreSeq(transcript['chrom'], transcript['strand'], exonEnds[ii], 
                                                    donor=True, verbose=verbose)])
        elif ii == transcript['exonCount'] - 1:
            acceptors = numpy.append(acceptors, [scoreSeq(transcript['chrom'], transcript['strand'], exonStarts[ii], 
                                                          donor=False, verbose=verbose)])
        else:
            donors = numpy.append(donors, [scoreSeq(transcript['chrom'], transcript['strand'], exonEnds[ii], 
                                                    donor=True, verbose=verbose)])
            acceptors = numpy.append(acceptors, [scoreSeq(transcript['chrom'], transcript['strand'], exonStarts[ii], 
                                                          donor=False, verbose=verbose)])
    return(donors, acceptors)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', "--outputFile", type=str)
    parser.add_argument('-t', "--transcripts", nargs='+', type=str)
    parser.add_argument('-v', "--verbose", type=bool, default=False)
    args = parser.parse_args()
    donors = numpy.array([])
    acceptors = numpy.array([])
    for this_transcript in args.transcripts:
        transcript = fetch_gene_coordinates(this_transcript)
        if args.verbose:
            print transcript
        transcript['exonStarts'] = re.sub(",(\s)*$", "", transcript['exonStarts'])
        transcript['exonEnds'] = re.sub(",(\s)*$", "", transcript['exonEnds'])
        (donors, acceptors) = addDataForThisTranscript(donors, acceptors, transcript, args.verbose)
    results = {
        "donors": {
            "mean": donors.mean(),
            "std": donors.std()
            },
        "acceptors": {
            "mean": acceptors.mean(),
            "std": acceptors.std()
            }
        }
    if args.verbose:
        print "donors mean", results["donors"]["mean"], "std", results["donors"]["std"]
        print "acceptors mean", results["acceptors"]["mean"], "std", results["acceptors"]["std"]
    with open(args.outputFile, "w") as fp:
        json.dump(results, fp)
        


if __name__ == "__main__":
    main()


