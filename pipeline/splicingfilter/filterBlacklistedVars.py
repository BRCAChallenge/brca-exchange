#!/usr/bin/env python
import subprocess

import click
import pandas as pd
import csv
from sys import platform

# priors-relevant columns that we'll initialize to '-' in the output
priorsCols = [
    "applicablePrior",
    "applicableEnigmaClass",
    "proteinPrior",
    "refDonorPrior",
    "deNovoDonorPrior",
    "refRefDonorMES",
    "refRefDonorZ",
    "altRefDonorMES",
    "altRefDonorZ",
    "refRefDonorSeq",
    "altRefDonorSeq",
    "refDonorVarStart",
    "refDonorVarLength",
    "refDonorExonStart",
    "refDonorIntronStart",
    "refDeNovoDonorMES",
    "refDeNovoDonorZ",
    "altDeNovoDonorMES",
    "altDeNovoDonorZ",
    "refDeNovoDonorSeq",
    "altDeNovoDonorSeq",
    "deNovoDonorVarStart",
    "deNovoDonorVarLength",
    "deNovoDonorExonStart",
    "deNovoDonorIntronStart",
    "deNovoDonorGenomicSplicePos",
    "deNovoDonorTranscriptSplicePos",
    "closestDonorGenomicSplicePos",
    "closestDonorTranscriptSplicePos",
    "closestDonorRefMES",
    "closestDonorRefZ",
    "closestDonorRefSeq",
    "closestDonorAltMES",
    "closestDonorAltZ",
    "closestDonorAltSeq",
    "closestDonorExonStart",
    "closestDonorIntronStart",
    "deNovoDonorAltGreaterRefFlag",
    "deNovoDonorAltGreaterClosestRefFlag",
    "deNovoDonorAltGreaterClosestAltFlag",
    "deNovoDonorFrameshiftFlag",
    "refAccPrior",
    "deNovoAccPrior",
    "refRefAccMES",
    "refRefAccZ",
    "altRefAccMES",
    "altRefAccZ",
    "refRefAccSeq",
    "altRefAccSeq",
    "refAccVarStart",
    "refAccVarLength",
    "refAccExonStart",
    "refAccIntronStart",
    "refDeNovoAccMES",
    "refDeNovoAccZ",
    "altDeNovoAccMES",
    "altDeNovoAccZ",
    "refDeNovoAccSeq",
    "altDeNovoAccSeq",
    "deNovoAccVarStart",
    "deNovoAccVarLength",
    "deNovoAccExonStart",
    "deNovoAccIntronStart",
    "deNovoAccGenomicSplicePos",
    "deNovoAccTranscriptSplicePos",
    "closestAccGenomicSplicePos",
    "closestAccTranscriptSplicePos",
    "closestAccRefMES",
    "closestAccRefZ",
    "closestAccRefSeq",
    "closestAccAltMES",
    "closestAccAltZ",
    "closestAccAltSeq",
    "closestAccExonStart",
    "closestAccIntronStart",
    "deNovoAccAltGreaterRefFlag",
    "deNovoAccAltGreaterClosestRefFlag",
    "deNovoAccAltGreaterClosestAltFlag",
    "deNovoAccFrameshiftFlag",
    "spliceSite",
    "spliceRescue",
    "spliceFlag",
    "frameshiftFlag",
    "inExonicPortionFlag",
    "CIDomainInRegionFlag",
    "isDivisibleFlag",
    "lowMESFlag",
    "varConsequences"
]


def filter_vars(variants, blacklisted_vars, output):
    """
    For every variant in the input, if the variant's RefSeq + HGVS_cDNA columns match an entry in
    blacklist, the values in that row for columns listed in priorsCols are replaced with '-'.
    :param variants: a file-like object containing variant information in TSV format
    :param blacklisted_vars: a file-like object containing a list of variant identifiers, like 'NM_000059.3:c.1909+1G>A'
    :param output: file-like object to which to write the resulting output
    :return: if output is None, the filtered output; otherwise, None
    """
    blacklist = set(x.strip() for x in blacklisted_vars.readlines())
    var_table = pd.read_csv(variants, sep='\t', na_filter=False)

    # for the blacklisted rows, initialize priors-relevant columns to '-'
    var_table.loc[
        (var_table['Reference_Sequence'] + ':' + var_table['HGVS_cDNA']).isin(blacklist),
        priorsCols
    ] = '-'
    return var_table.to_csv(path_or_buf=output, sep='\t', index=False, float_format='%s')


@click.group()
@click.option("--blacklisted_vars", type=click.File("r"), default="blacklisted_vars.txt",
              help="File of HGVS cDNA coordinates of variants which should return blank priors data")
@click.option("--output", type=click.File("w"), default="-",
              help="File to write filtered variants, defaults to stdout")
@click.pass_context
def cli(ctx, blacklisted_vars, output):
    ctx.obj = {
        "blacklisted_vars": blacklisted_vars,
        "output": output
    }

@cli.command("filter")
@click.argument("variants", type=click.File("r"))
@click.pass_context
def filter_vars_cmd(ctx, variants):
    filter_vars(variants, ctx.obj['blacklisted_vars'], ctx.obj['output'])


if __name__ == "__main__":
    cli()
