#!/usr/bin/env python

"""
calcVarPriors

Parses a tsv file (default built.tsv) containing variant information and for each variant in file
calculates either the prior probability of pathogenicity or a prior ENGIMA classification based on variant type and variant location
"""

import csv
from itertools import chain

import subprocess
import click
import traceback
import multiprocessing
import pytest
import pyhgvs.utils as pyhgvs_utils

from calc_priors.constants import BRCA1_RefSeq, BRCA2_RefSeq
from calc_priors.extract import getVarType
from calc_priors.compute import getVarLocation
from calc_priors.verify import getVarStrand
from calc_priors.priors import getPriorProbAfterGreyZoneSNS, getPriorProbSpliceDonorSNS, getPriorProbSpliceAcceptorSNS, \
    getPriorProbInGreyZoneSNS, getPriorProbInExonSNS, getPriorProbOutsideTranscriptBoundsSNS, getPriorProbInIntronSNS, \
    getPriorProbInUTRSNS

from calc_priors.utils import Benchmark, approximate_compare_tsv, MismatchException


def getVarData(variant, boundaries, variantData, genome, transcript):
    """
    Given variant, boundaries (either "priors" or "enigma') and list of dictionaries with variant data
    Genome is a SequenceFileDB for genome and transcript is a pyhgvs transcript object)
       both genome and transcript are necessary to convert from genomic to transcript coordinates
    Checks that variant is a single nucleotide substitution
    Determines prior prob dictionary based on variant location
    Return dictionary containing all values for all new prior prob fields
    """
    varLoc = getVarLocation(variant, boundaries)
    varType = getVarType(variant)
    blankDict = {"applicablePrior": "-",
                 "applicableEnigmaClass": "-",
                 "proteinPrior": "-",
                 "refDonorPrior": "-",
                 "deNovoDonorPrior": "-",
                 "refRefDonorMES": "-",
                 "refRefDonorZ": "-",
                 "altRefDonorMES": "-",
                 "altRefDonorZ": "-",
                 "refRefDonorSeq": "-",
                 "altRefDonorSeq": "-",
                 "refDonorVarStart": "-",
                 "refDonorVarLength": "-",
                 "refDonorExonStart": "-",
                 "refDonorIntronStart": "-",
                 "refDeNovoDonorMES": "-",
                 "refDeNovoDonorZ": "-",
                 "altDeNovoDonorMES": "-",
                 "altDeNovoDonorZ": "-",
                 "refDeNovoDonorSeq": "-",
                 "altDeNovoDonorSeq": "-",
                 "deNovoDonorVarStart": "-",
                 "deNovoDonorVarLength": "-",
                 "deNovoDonorExonStart": "-",
                 "deNovoDonorIntronStart": "-",
                 "deNovoDonorGenomicSplicePos": "-",
                 "deNovoDonorTranscriptSplicePos": "-",
                 "closestDonorGenomicSplicePos": "-",
                 "closestDonorTranscriptSplicePos": "-",
                 "closestDonorRefMES": "-",
                 "closestDonorRefZ": "-",
                 "closestDonorRefSeq": "-",
                 "closestDonorAltMES": "-",
                 "closestDonorAltZ": "-",
                 "closestDonorAltSeq": "-",
                 "closestDonorExonStart": "-",
                 "closestDonorIntronStart": "-",
                 "deNovoDonorAltGreaterRefFlag": "-",
                 "deNovoDonorAltGreaterClosestRefFlag": "-",
                 "deNovoDonorAltGreaterClosestAltFlag": "-",
                 "deNovoDonorFrameshiftFlag": "-",
                 "refAccPrior": "-",
                 "deNovoAccPrior": "-",
                 "refRefAccMES": "-",
                 "refRefAccZ": "-",
                 "altRefAccMES": "-",
                 "altRefAccZ": "-",
                 "refRefAccSeq": "-",
                 "altRefAccSeq": "-",
                 "refAccVarStart": "-",
                 "refAccVarLength": "-",
                 "refAccExonStart": "-",
                 "refAccIntronStart": "-",
                 "refDeNovoAccMES": "-",
                 "refDeNovoAccZ": "-",
                 "altDeNovoAccMES": "-",
                 "altDeNovoAccZ": "-",
                 "refDeNovoAccSeq": "-",
                 "altDeNovoAccSeq": "-",
                 "deNovoAccVarStart": "-",
                 "deNovoAccVarLength": "-",
                 "deNovoAccExonStart": "-",
                 "deNovoAccIntronStart": "-",
                 "deNovoAccGenomicSplicePos": "-",
                 "deNovoAccTranscriptSplicePos": "-",
                 "closestAccGenomicSplicePos": "-",
                 "closestAccTranscriptSplicePos": "-",
                 "closestAccRefMES": "-",
                 "closestAccRefZ": "-",
                 "closestAccRefSeq": "-",
                 "closestAccAltMES": "-",
                 "closestAccAltZ": "-",
                 "closestAccAltSeq": "-",
                 "closestAccExonStart": "-",
                 "closestAccIntronStart": "-",
                 "deNovoAccAltGreaterRefFlag": "-",
                 "deNovoAccAltGreaterClosestRefFlag": "-",
                 "deNovoAccAltGreaterClosestAltFlag": "-",
                 "deNovoAccFrameshiftFlag": "-",
                 "spliceSite": "-",
                 "spliceRescue": "-",
                 "spliceFlag": "-",
                 "frameshiftFlag": "-",
                 "inExonicPortionFlag": "-",
                 "CIDomainInRegionFlag": "-",
                 "isDivisibleFlag": "-",
                 "lowMESFlag": "-"}
    if varType == "substitution":
        # functions only work for variants with cannonical nucleotides (ACTG)
        if variant["Ref"] in ["A", "C", "G", "T"] and variant["Alt"] in ["A", "C", "G", "T"]:
            if varLoc == "outside_transcript_boundaries_variant":
                varData = getPriorProbOutsideTranscriptBoundsSNS(variant, boundaries)
            elif varLoc == "CI_splice_donor_variant" or varLoc == "splice_donor_variant":
                varData = getPriorProbSpliceDonorSNS(variant, boundaries, variantData, genome, transcript)
            elif varLoc == "CI_splice_acceptor_variant" or varLoc == "splice_acceptor_variant":
                varData = getPriorProbSpliceAcceptorSNS(variant, boundaries, variantData, genome, transcript)
            elif varLoc == "CI_domain_variant" or varLoc == "exon_variant":
                varData = getPriorProbInExonSNS(variant, boundaries, variantData, genome, transcript)
            elif varLoc == "grey_zone_variant":
                varData = getPriorProbInGreyZoneSNS(variant, boundaries, variantData)
            elif varLoc == "after_grey_zone_variant":
                varData = getPriorProbAfterGreyZoneSNS(variant, boundaries)
            elif varLoc == "UTR_variant":
                varData = getPriorProbInUTRSNS(variant, boundaries, genome, transcript)
            elif varLoc == "intron_variant":
                varData = getPriorProbInIntronSNS(variant, boundaries, genome, transcript)
            else:
                varData = blankDict.copy()
        else:
            varData = blankDict.copy()
    else:
        # to account for any non SNS variants
        varData = blankDict.copy()

    varData["varType"] = varType
    varData["varLoc"] = varLoc
    return varData


def addVarDataToRow(varData, inputRow):
    """
    Given data about a particular variant and a row from input file,
    Returns row with appended data
    """
    for key in varData.keys():
        inputRow[key] = varData[key]
    return inputRow


def getVarDict(variant, boundaries):
    """
    Given input data, returns a dictionary containing information for each variant in input
    Dictionary key is variant HGVS_cDNA and value is a dictionary containing variant gene, variant chromosome,
    variant strand, variant genomic coordinate, variant type, and variant location
    """
    varStrand = getVarStrand(variant)
    varType = getVarType(variant)
    varLoc = getVarLocation(variant, boundaries)

    varDict = {"varGene": variant["Gene_Symbol"],
               "varChrom": variant["Chr"],
               "varStrand": varStrand,
               "varGenCoordinate": variant["Pos"],
               "varType": varType,
               "varLoc": varLoc,
               "varHGVScDNA": variant["HGVS_cDNA"]}

    return varDict


brca1Transcript = None
brca2Transcript = None
useOldFile = True

if useOldFile:
    with open("mod_res_dn_brca20160525.txt", "r") as combinedfile:
        variantData = dict(
            ((x['gene'], x['nthgvs']), x) for x in csv.DictReader(combinedfile, delimiter="\t")
        )
else:
    with open("references/HCI_AllPriorsReport_BRCA1_V20.txt", "r") as brca1file, \
         open("references/HCI_AllPriorsReport_BRCA2_V20.txt", "r") as brca2file:
        variantData = dict(
            ((x['gene'], x['nthgvs']), x) for x in chain(
                csv.DictReader(brca1file, delimiter="\t"),
                csv.DictReader(brca2file, delimiter="\t")
            )
        )


def calc_one(variant):
    global brca1Transcript, brca2Transcript, variantData

    try:
        # variantData = csv.DictReader(open("mod_res_dn_brca20160525.txt", "r"), delimiter="\t")
        if variant["Gene_Symbol"] == "BRCA1":
            varData = getVarData(variant, "enigma", variantData, None, brca1Transcript)
        elif variant["Gene_Symbol"] == "BRCA2":
            varData = getVarData(variant, "enigma", variantData, None, brca2Transcript)
        click.echo("{}:{}".format(variant["HGVS_cDNA"], varData["varLoc"]), err=True)
        return addVarDataToRow(varData, variant)
    except Exception as e:
        traceback.print_exc()
        print('')
        raise e


def calc_all(variants, priors, genome, transcripts, processes):
    global brca1Transcript, brca2Transcript

    inputData = csv.DictReader(variants, delimiter="\t")
    fieldnames = inputData.fieldnames
    newHeaders = open("headers.tsv", "r").read().split()
    for header in newHeaders:
        fieldnames.append(header)
    outputData = csv.DictWriter(priors, delimiter="\t", lineterminator="\n", fieldnames=fieldnames)
    outputData.writerow(dict((fn, fn) for fn in inputData.fieldnames))

    # read RefSeq transcripts
    transcripts = pyhgvs_utils.read_transcripts(transcripts)

    brca1Transcript = transcripts.get(BRCA1_RefSeq)
    brca2Transcript = transcripts.get(BRCA2_RefSeq)

    if processes > 1:
        # Create a pool of processes and calculate in parallel
        click.echo("Processing using {} processes".format(processes), err=True)
        pool = multiprocessing.Pool(processes)
        try:
            # Normal map has a bug if there is no timout that prevents Keyboard interrupts:
            # https://stackoverflow.com/questions/1408356/keyboard-interrupts-with-pythons-multiprocessing-pool/1408476#1408476

            # calc_one_partial = functools.partial(calc_one, brca1=brca1Transcript, brca2=brca2Transcript)
            calculatedVariants = pool.map_async(calc_one, list(inputData)).get(99999999)

            # Sort output as the order of p.map is not deterministic
            outputData.writerows(sorted(
                calculatedVariants,
                key=lambda d: "{0}:g.{1}:{2}>{3}".format(d["Chr"], d["Pos"], d["Ref"], d["Alt"])))
        except KeyboardInterrupt:
            pool.terminate()
    else:
        outputData.writerows(sorted(
            map(calc_one, inputData),
            key=lambda d: "{0}:g.{1}:{2}>{3}".format(d["Chr"], d["Pos"], d["Ref"], d["Alt"])
        ))


def run(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    while True:
        line = process.stdout.readline().rstrip()
        if line:
            click.echo(line, err=True)
        else:
            break


@click.group()
@click.option("--genome", type=click.Path(exists=False), default="/references/hg38.fa",
              help="Fasta file containing hg38 reference genome")
@click.option("--transcripts", type=click.File("r"), default="refseq_annotation.hg38.gp",
              help="RefSeq annotation hg38-based genepred file")
@click.option("--processes", type=int, default=8,
              help="Number of processes to use")
@click.pass_context
def cli(ctx, genome, transcripts, processes):
    ctx.obj = {"genome": genome, "transcripts": transcripts, "processes": processes}


@cli.command()
def references():
    run("./references.sh")


@cli.command(help="Run only the unit tests")
@click.pass_context
def unittest(ctx):
    with Benchmark("ran unit tests"):
        pytest.main(["-p", "no:cacheprovider", "-x", "."])


@cli.command(help="Run self test")
@click.argument("length", type=click.Choice(["short", "long", "concerning"]))
@click.pass_context
def test(ctx, length):
    with Benchmark("running %s, time" % length):
        pytest.main(["-p", "no:cacheprovider", "-x", "."])
        calc_all(click.open_file("tests/variants_%s.tsv" % length, mode="r"),
                 click.open_file("/tmp/priors_%s.tsv" % length, mode="w"),
                 ctx.obj["genome"], ctx.obj["transcripts"], ctx.obj["processes"])

    # FIXME: we unfortunately can't use the md5sum test until we've replaced the sums
    # 1) we changed the protein_priors resource file
    # 2) we may switch to maxentpy, which produces slightly different values
    # for now, we defer to a slower approximate matching technique
    run("md5sum -c ./tests/md5/priors_%s.md5" % length)

    # now, run the test
    with Benchmark("test running %s, time" % length):
        infile = "tests/priors_%s.tsv" % length
        outfile = "/tmp/priors_%s.tsv" % length
        try:
            approximate_compare_tsv(infile, outfile)
        except MismatchException as ex:
            print("Mismatch between %s and %s; reason: %s" % (infile, outfile, str(ex)))


@cli.command()
@click.argument("variants", type=click.File("r"))
@click.argument("priors", type=click.File("w"))
@click.pass_context
def calc(ctx, variants, priors):
    calc_all(variants, priors,
             ctx.obj["genome"], ctx.obj["transcripts"], ctx.obj["processes"])


if __name__ == "__main__":
    cli()
