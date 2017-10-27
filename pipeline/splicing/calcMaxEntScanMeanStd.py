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

db = None

def fetch_gene_coordinates(transcript_name):
    """Query the indicated transcripts from the genome browser hg38 refseq table"""
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

def runMaxEntScan(sequence, donor=False):
    """Run maxEntScan on the indicated sequence.  Run score5.pl on candidate donor sequences, 
       score3.pl on candidate acceptor sequences.  Return the score"""
    (fd, tmpfile) = tempfile.mkstemp()
    fp = os.fdopen(fd, "w")
    fp.write(sequence)
    fp.close()
    if donor:
        pipe = subprocess.Popen(["perl","score5.pl", tmpfile], stdout=subprocess.PIPE)
    else:
        pipe = subprocess.Popen(["perl","score3.pl", tmpfile], stdout=subprocess.PIPE)
    result = pipe.stdout.read()
    entScore = re.findall("[+-]?\d+(?:\.\d+)?", str(result))
    os.remove(tmpfile)
    return(float(entScore[0]))


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


