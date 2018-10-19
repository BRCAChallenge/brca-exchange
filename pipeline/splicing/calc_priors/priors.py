from calc_priors.constants import PATHOGENIC_PROBABILITY, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH, \
    STD_ACC_INTRONIC_LENGTH, STD_ACC_EXONIC_LENGTH, LOW_MES_CUTOFF, HIGH_MES_CUTOFF, CAPPED_PROBABILITY, \
    LOW_SPLICING_PROBABILITY, REF_DONOR_REFZ_CUTOFF, MIN_REF_ALT_ZDIFF, HIGH_SPLICING_PROBABILITY, \
    MODERATE_SPLICING_PROBABILITY, REF_DONOR_HIGH_CUTOFF, REF_DONOR_LOW_CUTOFF, REF_ACC_REFZ_CUTOFF, \
    REF_ACC_HIGH_CUTOFF, REF_ACC_LOW_CUTOFF, LOW_PROBABILITY, DE_NOVO_DONOR_LOW_CUTOFF, DE_NOVO_DONOR_HIGH_CUTOFF, \
    MODERATE_DENOVO_PROBABILITY, HIGH_PROBABILITY, STD_DE_NOVO_OFFSET, STD_ACC_SIZE


from calc_priors import verify
from calc_priors import extract
from calc_priors import compute


def getPriorProbSpliceRescueNonsenseSNS(variant, boundaries, deNovoDonorInRefAcc=False):
    """
    Given a variant, determines if there is a possibility of splice rescue
    deNovoDonorInRefAcc argument = True  if looking for deNovoDonor in ref acceptor site, False otherwise
    If there is a possibility of splice rescue, flags variant for further analysis
    Else assigns prior probability of pathogenecity and predicted qualitative ENIGMA class
    Extra flags are included to provide more information about why splicing rescue does or does not occur
      frameshiftFlag: equals 1 if variant causes a frameshift, 0 if not
      inExonicPortionFlag: equals 1 if variant IS in exonic portion of highest scoring window, 0 if not
      CIDomainRegionFlag: equals 1 if truncating region (between new donor and next reference splice acceptor)
         includes a clinicall import domain, 0 if not
      isDivisibleFlag: equals 1 if distance between de novo donor and wild-type donor splice site is NOT divisible by 3,
         0 if not
    """
    if compute.varInIneligibleDeNovoExon(variant, donor=True):
        return {"priorProb": PATHOGENIC_PROBABILITY,
                "enigmaClass": "class_5",
                "spliceRescue": 0,
                "spliceFlag": 0,
                "frameshiftFlag": "N/A",
                "inExonicPortionFlag": "N/A",
                "CIDomainInRegionFlag": "N/A",
                "isDivisibleFlag": "N/A",
                "lowMESFlag": "N/A",
                "varConsequences": "N/A"}
    # checks that variant causes a premature stop codon in an exon
    varCons = extract.getVarConsequences(variant)
    if "stop_gained" in varCons and compute.varInExon(variant):
        spliceFlag = 0
        spliceRescue = 0
        frameshiftFlag = "-"
        inExonicPortionFlag = "-"
        CIDomainInRegionFlag = "-"
        isDivisibleFlag = "-"
        lowMESFlag = "-"
        # if variant is in specified exonic portion of highest scoring sliding window, no splice rescue
        if compute.varInExonicPortion(variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH, donor=True,
                              deNovoDonorInRefAcc=deNovoDonorInRefAcc):
            priorProb = PATHOGENIC_PROBABILITY
            inExonicPortionFlag = 1
            spliceRescue = 0
        else:
            inFrame = compute.isSplicingWindowInFrame(variant, STD_EXONIC_PORTION, STD_ACC_INTRONIC_LENGTH,
                                              deNovoDonorInRefAcc=deNovoDonorInRefAcc, donor=True)
            # if variant causes a frameshift, no splice rescue
            if not inFrame:
                priorProb = PATHOGENIC_PROBABILITY
                frameshiftFlag = 1
                inExonicPortionFlag = 0
                spliceRescue = 0
            else:
                varGenPos = int(variant["Pos"])
                varStrand = verify.getVarStrand(variant)
                varExonNum = compute.getVarExonNumberSNS(variant)
                # varExonNum returns a string in the format "exonN"
                # nextExonNum parses out N from varExonNum and adds 1 to get next exon number key "exonN+1"
                # use [4:] to remove "exon" from "exonN" so can add 1 to N to get N+1
                nextExonNum = "exon" + str(int(varExonNum[4:]) + 1)
                # skips to exon 5 because exon 4 does not exist in BRCA1 refseq transcript
                if variant["Gene_Symbol"] == "BRCA1" and nextExonNum == "exon4":
                    nextExonNum = "exon5"
                refSpliceAccBounds = extract.getSpliceAcceptorBoundaries(variant, STD_ACC_INTRONIC_LENGTH,
                                                                 STD_ACC_EXONIC_LENGTH)
                varWindowPos = compute.getVarWindowPosition(variant, donor=True, deNovoDonorInRefAcc=deNovoDonorInRefAcc)
                inExonicPortion = compute.varInExonicPortion(variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH, donor=True,
                                                     deNovoDonorInRefAcc=deNovoDonorInRefAcc)
                # gets region from new splice position to next splice acceptor
                regionStart = extract.getNewSplicePosition(varGenPos, varStrand, varWindowPos, inExonicPortion,
                                                   STD_EXONIC_PORTION,
                                                   STD_ACC_INTRONIC_LENGTH, donor=True)
                regionEnd = refSpliceAccBounds[nextExonNum]["acceptorStart"]
                CIDomainInRegion = verify.isCIDomainInRegion(regionStart, regionEnd, boundaries, variant["Gene_Symbol"])
                isDivisible = compute.isDeNovoWildTypeSplicePosDistanceDivisibleByThree(variant, STD_EXONIC_PORTION,
                                                                                STD_ACC_INTRONIC_LENGTH,
                                                                                deNovoDonorInRefAcc=deNovoDonorInRefAcc,
                                                                                donor=True)
                # if truncated region includes a clinically important domain or causes a frameshift
                if CIDomainInRegion:
                    priorProb = PATHOGENIC_PROBABILITY
                    spliceRescue = 0
                    CIDomainInRegionFlag = 1
                    inExonicPortionFlag = 0
                    frameshiftFlag = 0
                elif not isDivisible:
                    priorProb = PATHOGENIC_PROBABILITY
                    spliceRescue = 0
                    isDivisibleFlag = 1
                    frameshiftFlag = 1
                    inExonicPortionFlag = 0
                    CIDomainInRegionFlag = 0
                else:
                    # possibility of splice rescue, check MES scores
                    frameshiftFlag = 0
                    inExonicPortionFlag = 0
                    CIDomainInRegionFlag = 0
                    isDivisibleFlag = 0
                    deNovoSpliceData = compute.getMaxMaxEntScanScoreSlidingWindowSNS(variant, STD_EXONIC_PORTION,
                                                                             STD_DE_NOVO_LENGTH,
                                                                             donor=True,
                                                                             deNovoDonorInRefAcc=deNovoDonorInRefAcc)
                    altMES = deNovoSpliceData["altMaxEntScanScore"]
                    refZScore = deNovoSpliceData["refZScore"]
                    altZScore = deNovoSpliceData["altZScore"]
                    if altZScore <= refZScore:
                        # variant creates a weaker splice site than was there previously, no change in splicing
                        priorProb = PATHOGENIC_PROBABILITY
                        spliceRescue = 0
                        lowMESFlag = 1
                    else:
                        # variant creates a stronger splice site than reference sequence
                        if altMES < LOW_MES_CUTOFF:
                            # still a weak splice site, so no change in splicing
                            priorProb = PATHOGENIC_PROBABILITY
                            spliceRescue = 0
                            lowMESFlag = 1
                        elif (altMES >= LOW_MES_CUTOFF) and (altMES <= HIGH_MES_CUTOFF):
                            # still a weak splice site, but possibility of splice rescue
                            deNovoOffset = 0
                            closestRefData = compute.getClosestSpliceSiteScores(variant, deNovoOffset, donor=True, deNovo=False,
                                                                        deNovoDonorInRefAcc=deNovoDonorInRefAcc)
                            closestZScore = closestRefData["zScore"]
                            if compute.varInSpliceRegion(variant, donor=True, deNovo=False):
                                closestAltData = getPriorProbRefSpliceDonorSNS(variant, boundaries)
                                closestZScore = closestAltData["altZScore"]
                            if altZScore > closestZScore:
                                # splice site created by variant is stronger than subsequent (closest) wild-type donor
                                priorProb = "N/A"
                                enigmaClass = "N/A"
                                spliceRescue = 1
                                spliceFlag = 1
                                lowMESFlag = 0
                            else:
                                # splice site created by variant is weaker than subsequent (closest) wild-type donor
                                priorProb = PATHOGENIC_PROBABILITY
                                spliceRescue = 0
                                lowMESFlag = 1
                        else:
                            # altMES > HIGH_MES_CUTOFF, strong splice site with higher possibility of splicing rescue
                            priorProb = "N/A"
                            enigmaClass = "N/A"
                            spliceRescue = 1
                            spliceFlag = 1
                            lowMESFlag = 0

        # capped nonsense probability due to special cases of in-frame exon skipping
        # for BRCA1 exons 9/10 and BRCA2 exon 12 (only in splice donor/acceptor region)
        if spliceFlag == 0:
            # splice rescue does not occur
            varGene = variant["Gene_Symbol"]
            if compute.varInExon(variant):
                exonName = compute.getVarExonNumberSNS(variant)
                if varGene == "BRCA1" and (exonName == "exon9" or exonName == "exon10"):
                    priorProb = CAPPED_PROBABILITY
                if varGene == "BRCA2" and exonName == "exon12":
                    inSpliceDonor = compute.varInSpliceRegion(variant, donor=True, deNovo=False)
                    inSpliceAcceptor = compute.varInSpliceRegion(variant, donor=False, deNovo=False)
                    if inSpliceDonor or inSpliceAcceptor:
                        priorProb = CAPPED_PROBABILITY

        if priorProb == "N/A":
            pass
        else:
            enigmaClass = extract.getEnigmaClass(priorProb)

        return {"priorProb": priorProb,
                "enigmaClass": enigmaClass,
                "spliceRescue": spliceRescue,
                "spliceFlag": spliceFlag,
                "frameshiftFlag": frameshiftFlag,
                "inExonicPortionFlag": inExonicPortionFlag,
                "CIDomainInRegionFlag": CIDomainInRegionFlag,
                "isDivisibleFlag": isDivisibleFlag,
                "lowMESFlag": lowMESFlag,
                "varConsequences": ",".join(varCons)}


def getPriorProbRefSpliceDonorSNS(variant, boundaries):
    """
    Given a variant and location boundaries (either PRIORS or enigma)
    Checks that variant is in a splice donor site and is a single nucleotide substitution
    Returns a dictionary containing:
     prior probability of pathogenecity, predicted qualitative enigma class, ref and alt sequences, and ref and alt MES and zscores
     also contains variables needed to format ref and alt seq:
       - varStart: index (based on 0 indexing) where variant is located in altSeq
       - varLength: length of variant in altSeq (this function is only for SNS variants so varLength=1)
       - exonStart: denotes index (0-indexed) where exonic portion of sequence starts, always equal to 0 because splice donor is at end of exon
       - intronStart: denotes index (0-indexed) where intronic portion of sequence starts, equal to exonic portion size
     also has spliceSite variable = 1, so variant is marked as in a reference splice site
    """
    varType = extract.getVarType(variant)
    varLoc = compute.getVarLocation(variant, boundaries)
    if varType == "substitution" and (varLoc == "splice_donor_variant" or varLoc == "CI_splice_donor_variant"):
        # to get region boundaries to get ref and alt seq
        spliceDonorBounds = compute.getVarSpliceRegionBounds(variant, donor=True, deNovo=False)
        refAltSeqs = extract.getRefAltSeqs(variant, spliceDonorBounds["donorStart"], spliceDonorBounds["donorEnd"])
        scores = extract.getRefAltScores(refAltSeqs["refSeq"], refAltSeqs["altSeq"], donor=True)
        refMaxEntScanScore = scores["refScores"]["maxEntScanScore"]
        refZScore = scores["refScores"]["zScore"]
        altMaxEntScanScore = scores["altScores"]["maxEntScanScore"]
        altZScore = scores["altScores"]["zScore"]
        if altMaxEntScanScore >= refMaxEntScanScore:
            priorProb = LOW_SPLICING_PROBABILITY
        elif (refZScore < REF_DONOR_REFZ_CUTOFF) and ((refZScore - altZScore) > MIN_REF_ALT_ZDIFF):
            priorProb = HIGH_SPLICING_PROBABILITY
        elif (refZScore < REF_DONOR_REFZ_CUTOFF) and ((refZScore - altZScore) < MIN_REF_ALT_ZDIFF):
            priorProb = MODERATE_SPLICING_PROBABILITY
        else:
            if altZScore > REF_DONOR_HIGH_CUTOFF:
                priorProb = LOW_SPLICING_PROBABILITY
            elif altZScore <= REF_DONOR_HIGH_CUTOFF and altZScore >= REF_DONOR_LOW_CUTOFF:
                priorProb = MODERATE_SPLICING_PROBABILITY
            else:
                priorProb = HIGH_SPLICING_PROBABILITY

        # capped splicing probability due to special cases of in-frame exon skipping
        varGene = variant["Gene_Symbol"]
        exonName = spliceDonorBounds["exonName"]
        if varGene == "BRCA1" and (exonName == "exon9" or exonName == "exon10"):
            if priorProb == HIGH_SPLICING_PROBABILITY:
                priorProb = CAPPED_PROBABILITY
        if varGene == "BRCA2" and exonName == "exon12":
            if priorProb == HIGH_SPLICING_PROBABILITY:
                priorProb = CAPPED_PROBABILITY

        enigmaClass = extract.getEnigmaClass(priorProb)

        return {"priorProb": priorProb,
                "enigmaClass": enigmaClass,
                "refSeq": refAltSeqs["refSeq"].upper(),
                "altSeq": refAltSeqs["altSeq"].upper(),
                "varStart": extract.getVarSeqIndexSNS(refAltSeqs["refSeq"], refAltSeqs["altSeq"]),
                "varLength": 1,
                "exonStart": 0,
                "intronStart": STD_EXONIC_PORTION,
                "refMaxEntScanScore": refMaxEntScanScore,
                "altMaxEntScanScore": altMaxEntScanScore,
                "refZScore": refZScore,
                "altZScore": altZScore,
                "spliceSite": 1}


def getPriorProbRefSpliceAcceptorSNS(variant, boundaries):
    """
    Given a variant and location boundaries (either PRIORS or enigma)
    Checks that variant is in a splice acceptor site and is a single nucleotide substitution
    Returns a dictionary containing:
     prior probability of pathogenecity, predicted qualitative enigma class, ref and alt sequences, and ref and alt MES and zscores
     also contains variables needed to format ref and alt seq:
       - varStart: index (based on 0 indexing) where variant is located in altSeq
       - varLength: length of variant in altSeq (this function is only for SNS variants so varLength=1)
       - exonStart: denotes index (0-indexed) where exonic portion of sequence starts, always equal to 0 because splice donor is at end of exon
       - intronStart: denotes index (0-indexed) where intronic portion of sequence starts, equal to exonic portion size
     also has spliceSite variable = 1, so variant is marked as in a reference splice site
    """
    varType = extract.getVarType(variant)
    varLoc = compute.getVarLocation(variant, boundaries)
    if varType == "substitution" and (varLoc == "splice_acceptor_variant" or varLoc == "CI_splice_acceptor_variant"):
        # to get region boundaires to get ref and alt seq
        spliceAcceptorBounds = compute.getVarSpliceRegionBounds(variant, donor=False, deNovo=False)
        refAltSeqs = extract.getRefAltSeqs(variant, spliceAcceptorBounds["acceptorStart"], spliceAcceptorBounds["acceptorEnd"])
        scores = extract.getRefAltScores(refAltSeqs["refSeq"], refAltSeqs["altSeq"], donor=False)
        refMaxEntScanScore = scores["refScores"]["maxEntScanScore"]
        refZScore = scores["refScores"]["zScore"]
        altMaxEntScanScore = scores["altScores"]["maxEntScanScore"]
        altZScore = scores["altScores"]["zScore"]
        if altMaxEntScanScore >= refMaxEntScanScore:
            priorProb = LOW_SPLICING_PROBABILITY
        elif (refZScore < REF_ACC_REFZ_CUTOFF) and ((refZScore - altZScore) > MIN_REF_ALT_ZDIFF):
            priorProb = HIGH_SPLICING_PROBABILITY
        elif (refZScore < REF_ACC_REFZ_CUTOFF) and ((refZScore - altZScore) < MIN_REF_ALT_ZDIFF):
            priorProb = MODERATE_SPLICING_PROBABILITY
        else:
            if altZScore > REF_ACC_HIGH_CUTOFF:
                priorProb = LOW_SPLICING_PROBABILITY
            elif altZScore <= REF_ACC_HIGH_CUTOFF and altZScore >= REF_ACC_LOW_CUTOFF:
                priorProb = MODERATE_SPLICING_PROBABILITY
            else:
                priorProb = HIGH_SPLICING_PROBABILITY

        # capped splicing probability due to special cases of in-frame exon skipping
        varGene = variant["Gene_Symbol"]
        exonName = spliceAcceptorBounds["exonName"]
        if varGene == "BRCA1" and (exonName == "exon9" or exonName == "exon10"):
            if priorProb == HIGH_SPLICING_PROBABILITY:
                priorProb = CAPPED_PROBABILITY
        if varGene == "BRCA2" and exonName == "exon12":
            if priorProb == HIGH_SPLICING_PROBABILITY:
                priorProb = CAPPED_PROBABILITY

        enigmaClass = extract.getEnigmaClass(priorProb)

        return {"priorProb": priorProb,
                "enigmaClass": enigmaClass,
                "refSeq": refAltSeqs["refSeq"].upper(),
                "altSeq": refAltSeqs["altSeq"].upper(),
                "varStart": extract.getVarSeqIndexSNS(refAltSeqs["refSeq"], refAltSeqs["altSeq"]),
                "varLength": 1,
                "exonStart": len(refAltSeqs["refSeq"]) - STD_EXONIC_PORTION,
                "intronStart": 0,
                "refMaxEntScanScore": refMaxEntScanScore,
                "altMaxEntScanScore": altMaxEntScanScore,
                "refZScore": refZScore,
                "altZScore": altZScore,
                "spliceSite": 1}


def getPriorProbAfterGreyZoneSNS(variant, boundaries):
    """
    Given a variant and location boundaries (either PRIORS or enigma)
    Checks that variant is after the grey zone and is a single nucleotide substitution
    Checks that variant is either a missense or nonsense mutation
    Returns a dictionary containing prior probability of pathogenecity and predicted qualitative enigma class
    Dictionary contains other values that are set to either N/A, "-", or 0 because they are not relevant
    """
    varType = extract.getVarType(variant)
    varLoc = compute.getVarLocation(variant, boundaries)
    if varType == "substitution" and varLoc == "after_grey_zone_variant":
        varCons = extract.getVarConsequences(variant)
        if "stop_gained" in varCons or "missense_variant" in varCons or "synonymous_variant" in varCons:
            priorProb = "N/A"
            enigmaClass = "class_2"
        return {"applicablePrior": priorProb,
                "applicableEnigmaClass": enigmaClass,
                "proteinPrior": "N/A",
                "refDonorPrior": "N/A",
                "deNovoDonorPrior": "N/A",
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
                "closestDonorRefMES": "-",
                "closestDonorRefZ": "-",
                "closestDonorRefSeq": "-",
                "closestDonorAltMES": "-",
                "closestDonorAltZ": "-",
                "closestDonorAltSeq": "-",
                "closestDonorExonStart": "-",
                "closestDonorIntronStart": "-",
                "deNovoDonorAltGreaterRefFlag": "N/A",
                "deNovoDonorAltGreaterClosestRefFlag": "N/A",
                "deNovoDonorAltGreaterClosestAltFlag": "N/A",
                "deNovoDonorFrameshiftFlag": "N/A",
                "refAccPrior": "N/A",
                "deNovoAccPrior": "N/A",
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
                "closestAccRefMES": "-",
                "closestAccRefZ": "-",
                "closestAccRefSeq": "-",
                "closestAccAltMES": "-",
                "closestAccAltZ": "-",
                "closestAccAltSeq": "-",
                "closestAccExonStart": "-",
                "closestAccIntronStart": "-",
                "deNovoAccAltGreaterRefFlag": "N/A",
                "deNovoAccAltGreaterClosestRefFlag": "N/A",
                "deNovoAccAltGreaterClosestAltFlag": "N/A",
                "deNovoAccFrameshiftFlag": "N/A",
                "spliceSite": 0,
                "spliceRescue": "N/A",
                "spliceFlag": 0,
                "frameshiftFlag": "N/A",
                "inExonicPortionFlag": "N/A",
                "CIDomainInRegionFlag": "N/A",
                "isDivisibleFlag": "N/A",
                "lowMESFlag": "N/A",
                "varConsequences": ",".join(varCons)}

    assert False, "Should never reach this"


def getPriorProbDeNovoDonorSNS(variant, boundaries, exonicPortionSize, genome, transcript, deNovoDonorInRefAcc=False):
    """
    Given a variant, boundaries (either priors or enigma), and exonicPortionSize
      1. checks that variant is a single nucleotide substitution
      2. checks that variant is in an exon or is in a reference splice donor region
    Genome is a SequenceFileDB for genome and transcript is a pyhgvs transcript object)
       both genome and transcript are necessary to convert from genomic to transcript coordinates
    Returns a dictionary containing:
      prior probability of pathogenecity and predicted qualitative engima class
      deNovo donor MaxEntScan scores, zscores, and sequences for ref and alt
      deNovo donor genomic and transcript splice positions in format "g.N" and "c.N" respectively
      closest donor MaxEntScan score, zscore, and sequence
      closest donor genomic and transcript splice positions in format "g.N" and "c.N" respectively
      altGreaterRefFlag: equals 1 if altZScore > refZScore, 0 otherwise
      altGreaterClosestRefFlag: equals 1 if altZScore > closestRefZScore, 0 otherwise
      altGreaterClosestAltFlag: equals 1 if variant is in a native donor and altZScore > closestAltZScore, 0 if not,
                                N/A if not in a native donor
      frameshiftFlag: equals 1 if new splice site causes a frameshift
      variables necessary for formatting any sequences that are returned by this function
    deNovoDonorInRefAcc = False if NOT checking for de novo donor in reference splice acceptor site
    deNovoDonorInRefAcc = True if checking for de novo donors in reference splice acceptor site
    """
    if extract.getVarType(variant) == "substitution":
        if compute.varInSpliceRegion(variant, donor=True, deNovo=True):
            if compute.varInExon(variant) and compute.varInIneligibleDeNovoExon(variant, donor=True):
                return {"priorProb": "N/A",
                        "enigmaClass": "N/A",
                        "refMaxEntScanScore": "N/A",
                        "altMaxEntScanScore": "N/A",
                        "refZScore": "N/A",
                        "altZScore": "N/A",
                        "refSeq": "N/A",
                        "altSeq": "N/A",
                        "varStart": "N/A",
                        "varLength": "N/A",
                        "exonStart": "N/A",
                        "intronStart": "N/A",
                        "genomicSplicePos": "N/A",
                        "transcriptSplicePos": "N/A",
                        "closestGenomicSplicePos": "N/A",
                        "closestTranscriptSplicePos": "N/A",
                        "closestRefMaxEntScanScore": "N/A",
                        "closestRefZScore": "N/A",
                        "closestRefSeq": "N/A",
                        "closestAltMaxEntScanScore": "N/A",
                        "closestAltZScore": "N/A",
                        "closestAltSeq": "N/A",
                        "closestExonStart": "N/A",
                        "closestIntronStart": "N/A",
                        "frameshiftFlag": "N/A",
                        "altGreaterRefFlag": "N/A",
                        "altGreaterClosestRefFlag": "N/A",
                        "altGreaterClosestAltFlag": "N/A"}
            slidingWindowInfo = compute.getMaxMaxEntScanScoreSlidingWindowSNS(variant, exonicPortionSize, STD_DE_NOVO_LENGTH,
                                                                      donor=True,
                                                                      deNovoDonorInRefAcc=deNovoDonorInRefAcc)
            # +-1 because per HCI PRIORS website donor position is defined as being first nucleotide that is NOT included in spliced exon
            if verify.getVarStrand(variant) == "-":
                newGenomicSplicePos = extract.getNewSplicePosition(variant["Pos"], "-", slidingWindowInfo["varWindowPosition"],
                                                           slidingWindowInfo["inExonicPortion"], STD_EXONIC_PORTION,
                                                           STD_ACC_INTRONIC_LENGTH,
                                                           donor=True) - 1
            else:
                newGenomicSplicePos = extract.getNewSplicePosition(variant["Pos"], "+", slidingWindowInfo["varWindowPosition"],
                                                           slidingWindowInfo["inExonicPortion"], STD_EXONIC_PORTION,
                                                           STD_ACC_INTRONIC_LENGTH,
                                                           donor=True) + 1
            deNovoOffset = 0
            subDonorInfo = compute.getClosestSpliceSiteScores(variant, deNovoOffset, donor=True, deNovo=False,
                                                      deNovoDonorInRefAcc=deNovoDonorInRefAcc)
            subZScore = subDonorInfo["zScore"]
            refAltZScore = "N/A"
            refAltMES = "N/A"
            refAltSeq = "N/A"
            altGreaterClosestRefFlag = 0
            altGreaterClosestAltFlag = "N/A"
            if compute.varInSpliceRegion(variant, donor=True, deNovo=False):
                refDonorInfo = getPriorProbRefSpliceDonorSNS(variant, "enigma")
                refAltZScore = refDonorInfo["altZScore"]
                refAltMES = refDonorInfo["altMaxEntScanScore"]
                refAltSeq = refDonorInfo["altSeq"]
                altGreaterClosestAltFlag = 0
            frameshiftFlag = 0
            frameshiftStatus = compute.getDeNovoSpliceFrameshiftStatus(variant, donor=True,
                                                               deNovoDonorInRefAcc=deNovoDonorInRefAcc)
            if frameshiftStatus:
                frameshiftFlag = 1
            altZScore = slidingWindowInfo["altZScore"]
            refZScore = slidingWindowInfo["refZScore"]
            altGreaterRefFlag = 0
            if altZScore <= refZScore:
                priorProb = LOW_PROBABILITY
            else:
                altGreaterRefFlag = 1
                if altZScore < DE_NOVO_DONOR_LOW_CUTOFF:
                    priorProb = LOW_PROBABILITY
                elif altZScore >= DE_NOVO_DONOR_LOW_CUTOFF and altZScore < DE_NOVO_DONOR_HIGH_CUTOFF:
                    priorProb = MODERATE_DENOVO_PROBABILITY
                else:
                    priorProb = HIGH_PROBABILITY
            if (altZScore > subZScore and refAltZScore == "N/A") or (
                    altZScore > refAltZScore and refAltZScore != "N/A"):
                # promote prior prob by one step
                if priorProb == LOW_PROBABILITY:
                    priorProb = MODERATE_DENOVO_PROBABILITY
                elif priorProb == MODERATE_DENOVO_PROBABILITY:
                    priorProb = HIGH_PROBABILITY
                else:
                    priorProb = priorProb

            # capped splicing probability due to special cases of in-frame exon skipping
            varGene = variant["Gene_Symbol"]
            exonName = subDonorInfo["exonName"]
            if varGene == "BRCA1" and (exonName == "exon9" or exonName == "exon10"):
                if priorProb == HIGH_PROBABILITY:
                    priorProb = CAPPED_PROBABILITY

            if altZScore > subZScore:
                altGreaterClosestRefFlag = 1
            if altZScore > refAltZScore and refAltZScore != "N/A":
                altGreaterClosestAltFlag = 1

            if frameshiftFlag == 0 and priorProb != 0:
                frameshiftAndCIStatus = compute.getDeNovoFrameshiftAndCIStatus(variant, boundaries, donor=True,
                                                                       deNovoDonorInRefAcc=deNovoDonorInRefAcc)
                if frameshiftAndCIStatus:
                    priorProb = LOW_PROBABILITY

            # converts genomic splice position to transcript splice position
            newTranscriptSplicePos = verify.convertGenomicPosToTranscriptPos(newGenomicSplicePos, extract.getVarChrom(variant), genome,
                                                                      transcript)
            # converts closest genomic splice position to transcript splice position
            closestGenomicSplicePos = subDonorInfo["genomicSplicePos"]
            closestTranscriptSplicePos = verify.convertGenomicPosToTranscriptPos(closestGenomicSplicePos, extract.getVarChrom(variant),
                                                                          genome, transcript)

            return {"priorProb": priorProb,
                    "enigmaClass": extract.getEnigmaClass(priorProb),
                    "refMaxEntScanScore": slidingWindowInfo["refMaxEntScanScore"],
                    "altMaxEntScanScore": slidingWindowInfo["altMaxEntScanScore"],
                    "refZScore": refZScore,
                    "altZScore": altZScore,
                    "refSeq": slidingWindowInfo["refSeq"],
                    "altSeq": slidingWindowInfo["altSeq"],
                    "varStart": slidingWindowInfo["varStart"],
                    "varLength": slidingWindowInfo["varLength"],
                    "exonStart": 0,
                    "intronStart": STD_EXONIC_PORTION,
                    "genomicSplicePos": verify.formatSplicePosition(newGenomicSplicePos, transcript=False),
                    "transcriptSplicePos": verify.formatSplicePosition(newTranscriptSplicePos, transcript=True),
                    "closestGenomicSplicePos": verify.formatSplicePosition(closestGenomicSplicePos, transcript=False),
                    "closestTranscriptSplicePos": verify.formatSplicePosition(closestTranscriptSplicePos, transcript=True),
                    "closestRefMaxEntScanScore": subDonorInfo["maxEntScanScore"],
                    "closestRefZScore": subDonorInfo["zScore"],
                    "closestRefSeq": subDonorInfo["sequence"],
                    "closestAltMaxEntScanScore": refAltMES,
                    "closestAltZScore": refAltZScore,
                    "closestAltSeq": refAltSeq,
                    "closestExonStart": subDonorInfo["exonStart"],
                    "closestIntronStart": subDonorInfo["intronStart"],
                    "altGreaterRefFlag": altGreaterRefFlag,
                    "altGreaterClosestRefFlag": altGreaterClosestRefFlag,
                    "altGreaterClosestAltFlag": altGreaterClosestAltFlag,
                    "frameshiftFlag": frameshiftFlag}


def getPriorProbDeNovoAcceptorSNS(variant, exonicPortionSize, deNovoLength, genome, transcript):
    """
    Given a variant, exonic portion size, and de novo length:
      1. checks that variant is a single nucleotide substitution
      2. checks that variant is in de novo splice acceptor region
         de novo splice acceptor region defined by deNovoLength
    Genome is a SequenceFileDB for genome and transcript is a pyhgvs transcript object)
       both genome and transcript are necessary to convert from genomic to transcript coordinates
    Returns a dictionary containing:
      prior probability of pathogenecity and predicted qualitative engima class (both N/A)
      deNovo acceptor MaxEntScan scores, zscores, and sequence for ref and alt
      de novo acceptor genomic and transcript splice positions in format "g.N" and "c.N" respectively
      MaxEntScan score, zscore, and sequence for closest ref splice acceptor
      closest ref acceptor genomic and transcript splice positions in format "g.N" and "c.N" respectively
      frameshiftFlag: equals 1 if new splice site causes a frameshift
      altGreaterRefFlag: which equals 1 altZScore > refZScore, 0 otherwise
      altGreaterClosestRefFlag: which equals 1 if altZScore > closestRefZScore, 0 otherwise
      altGreaterClosestAltFlag: which equals 1 if variant in native acceptor and altZScore > closestAltZScore, if not,
                                N/A if not in a native acceptor
      variables necessary for formatting the returned sequences
    """
    if extract.getVarType(variant) == "substitution":
        if compute.varInSpliceRegion(variant, donor=False, deNovo=True):
            slidingWindowInfo = compute.getMaxMaxEntScanScoreSlidingWindowSNS(variant, exonicPortionSize, deNovoLength,
                                                                      donor=False)
            newGenomicSplicePos = extract.getNewSplicePosition(variant["Pos"], verify.getVarStrand(variant),
                                                       slidingWindowInfo["varWindowPosition"],
                                                       slidingWindowInfo["inExonicPortion"], STD_EXONIC_PORTION,
                                                       STD_ACC_INTRONIC_LENGTH,
                                                       donor=False)
            deNovoOffset = deNovoLength - exonicPortionSize
            closestAccInfo = compute.getClosestSpliceSiteScores(variant, STD_DE_NOVO_OFFSET, donor=False, deNovo=True,
                                                        deNovoDonorInRefAcc=False)
            refAltZScore = "N/A"
            refAltMES = "N/A"
            refAltSeq = "N/A"
            altGreaterClosestRefFlag = 0
            altGreaterClosestAltFlag = "N/A"
            if compute.varInSpliceRegion(variant, donor=False, deNovo=False):
                refAccInfo = getPriorProbRefSpliceAcceptorSNS(variant, "enigma")
                refAltZScore = refAccInfo["altZScore"]
                refAltMES = refAccInfo["altMaxEntScanScore"]
                refAltSeq = refAccInfo["altSeq"]
                altGreaterClosestAltFlag = 0
            frameshiftFlag = 0
            frameshiftStatus = compute.getDeNovoSpliceFrameshiftStatus(variant, donor=False, deNovoDonorInRefAcc=False)
            if frameshiftStatus:
                frameshiftFlag = 1
            altZScore = slidingWindowInfo["altZScore"]
            refZScore = slidingWindowInfo["refZScore"]
            closestZScore = closestAccInfo["zScore"]
            altGreaterRefFlag = 0
            if altZScore > refZScore:
                altGreaterRefFlag = 1
            if altZScore > closestZScore:
                altGreaterClosestRefFlag = 1
            if altZScore > refAltZScore and refAltZScore != "N/A":
                altGreaterClosestAltFlag = 1

            # converts genomic splice position to transcript splice position
            newTranscriptSplicePos = verify.convertGenomicPosToTranscriptPos(newGenomicSplicePos, extract.getVarChrom(variant), genome,
                                                                      transcript)
            # converts closest genomic splice position to transcript splice position
            closestGenomicSplicePos = closestAccInfo["genomicSplicePos"]
            closestTranscriptSplicePos = verify.convertGenomicPosToTranscriptPos(closestGenomicSplicePos, extract.getVarChrom(variant),
                                                                          genome, transcript)

            return {"priorProb": "N/A",
                    "enigmaClass": "N/A",
                    "refMaxEntScanScore": slidingWindowInfo["refMaxEntScanScore"],
                    "altMaxEntScanScore": slidingWindowInfo["altMaxEntScanScore"],
                    "refZScore": refZScore,
                    "altZScore": altZScore,
                    "refSeq": slidingWindowInfo["refSeq"],
                    "altSeq": slidingWindowInfo["altSeq"],
                    "varStart": slidingWindowInfo["varStart"],
                    "varLength": slidingWindowInfo["varLength"],
                    "exonStart": STD_ACC_SIZE - STD_EXONIC_PORTION,
                    "intronStart": 0,
                    "genomicSplicePos": verify.formatSplicePosition(newGenomicSplicePos, transcript=False),
                    "transcriptSplicePos": verify.formatSplicePosition(newTranscriptSplicePos, transcript=True),
                    "closestGenomicSplicePos": verify.formatSplicePosition(closestGenomicSplicePos, transcript=False),
                    "closestTranscriptSplicePos": verify.formatSplicePosition(closestTranscriptSplicePos, transcript=True),
                    "closestRefMaxEntScanScore": closestAccInfo["maxEntScanScore"],
                    "closestRefZScore": closestAccInfo["zScore"],
                    "closestRefSeq": closestAccInfo["sequence"],
                    "closestAltMaxEntScanScore": refAltMES,
                    "closestAltZScore": refAltZScore,
                    "closestAltSeq": refAltSeq,
                    "closestExonStart": closestAccInfo["exonStart"],
                    "closestIntronStart": closestAccInfo["intronStart"],
                    "altGreaterRefFlag": altGreaterRefFlag,
                    "altGreaterClosestRefFlag": altGreaterClosestRefFlag,
                    "altGreaterClosestAltFlag": altGreaterClosestAltFlag,
                    "frameshiftFlag": frameshiftFlag}


def getPriorProbSpliceDonorSNS(variant, boundaries, variantData, genome, transcript):
    """
    Given a variant, boundaries (either PRIORS or ENIGMA), and a list of dictionaries with variant data
    Genome is a SequenceFileDB for genome and transcript is a pyhgvs transcript object)
       both genome and transcript are necessary to convert from genomic to transcript coordinates
    Determines reference donor and de novo donor scores for variant
    If variant causes a nonsense mutation, determines if splice rescue occurs
    Returns dicitionary containing scores and sequences for ref and de novo splice donor/acceptor,
        closest splice site scores and sequences, and protein prior if variant in exon
        "N/A" if score, sequence or flag not applicable for variant
    Also contains other values:
        applicable prior, highest prior if variant has multiple priors
        ref prior, prior prob for reference splice donor
        de novo prior, prior prob for de novo donor sequence
        splice flag = 1, because variant in reference splice site
        deNovoDonorFrameshiftFlag = 1 if new splice site causes a frameshift
        deNovoDonorAltGreaterRefFlag = 1 if variant is possible de novo donor variant
        deNovoDonorAltGreaterClosestRefFlag = 1 if variant is possible de novo donor variant
        deNovoDonorAltGreaterClosestAltFlag = 1 if variant is possible de novo donor variant
        deNovoAccFrameshiftFlag = N/A, because not applicable for variants in ref donor sites
        deNovoAccAltGreaterRefFlag = N/A, because not applicable for variants in ref donor sites
        deNovoAccAltGreaterClosestRefFlag = N/A, because not applicable for variants in ref donor sites
        deNovoAccAltGreaterClosestAltFlag = N/A, because not applicable for variants in ref donor sites
        spliceRescue = 1 if splice rescue possible for nonsense variant, 0 otherwise, N/A if not nonsense variant
        spliceFlag = 1 if splice rescue is possible so variant can be flagged for further analysis, 0 otherwise
        frameshiftFlag = 1 if nonsense variant causes a frameshift mutation also, 0 if not, N/A if not nonsense variant
        inExonicPortionFlag = 1 if variant IS in exonic portion of highest scoring window and variant is a nonsense variant,
           0 if NOT in exonic portion, N/A if not nonsense variant
        CIDomainRegionFlag = 1 if truncating region (between new donor and next reference splice acceptor)
           includes a clinically import domain for nonsense variant, 0 if not, N/A if not nonsense variant
        isDivisibleFlag = 1 if distance between de novo donor and wild-type donor splice site is NOT divisible by 3 for nonsense variant,
           0 if it IS divisible, N/A if not a nonsense variant
    Dictionary also contains formatting variables for each listed sequence
    """
    if compute.varInSpliceRegion(variant, donor=True, deNovo=False) and extract.getVarType(variant) == "substitution":
        refSpliceInfo = getPriorProbRefSpliceDonorSNS(variant, boundaries)
        deNovoSpliceInfo = getPriorProbDeNovoDonorSNS(variant, boundaries, STD_EXONIC_PORTION, genome, transcript,
                                                      deNovoDonorInRefAcc=False)
        deNovoPrior = deNovoSpliceInfo["priorProb"]
        refPrior = refSpliceInfo["priorProb"]
        proteinPrior = "N/A"
        if compute.varInExon(variant):
            proteinInfo = getPriorProbProteinSNS(variant, variantData)
            proteinPrior = proteinInfo["priorProb"]
        if deNovoPrior != "N/A" and proteinPrior != "N/A":
            applicablePrior = max(deNovoPrior, refPrior, proteinPrior)
        elif deNovoPrior == "N/A" and proteinPrior != "N/A":
            applicablePrior = max(refPrior, proteinPrior)
        elif deNovoPrior != "N/A" and proteinPrior == "N/A":
            applicablePrior = max(deNovoPrior, refPrior)
        elif deNovoPrior == "N/A" and proteinPrior == "N/A":
            applicablePrior = refPrior

        spliceRescue = "N/A"
        spliceFlag = 0
        frameshiftFlag = "N/A"
        inExonicPortionFlag = "N/A"
        CIDomainInRegionFlag = "N/A"
        isDivisibleFlag = "N/A"
        lowMESFlag = "N/A"
        # to check for nonsense variants in exonic portion of splice donor site
        varCons = extract.getVarConsequences(variant)
        if compute.varInExon(variant) and "stop_gained" in varCons:
            nonsenseData = getPriorProbSpliceRescueNonsenseSNS(variant, boundaries)
            applicablePrior = nonsenseData["priorProb"]
            spliceRescue = nonsenseData["spliceRescue"]
            spliceFlag = nonsenseData["spliceFlag"]
            frameshiftFlag = nonsenseData["frameshiftFlag"]
            inExonicPortionFlag = nonsenseData["inExonicPortionFlag"]
            CIDomainInRegionFlag = nonsenseData["CIDomainInRegionFlag"]
            isDivisibleFlag = nonsenseData["isDivisibleFlag"]
            lowMESFlag = nonsenseData["lowMESFlag"]

        return {"applicablePrior": applicablePrior,
                "applicableEnigmaClass": extract.getEnigmaClass(applicablePrior),
                "proteinPrior": proteinPrior,
                "refDonorPrior": refPrior,
                "deNovoDonorPrior": deNovoPrior,
                "refRefDonorMES": refSpliceInfo["refMaxEntScanScore"],
                "refRefDonorZ": refSpliceInfo["refZScore"],
                "altRefDonorMES": refSpliceInfo["altMaxEntScanScore"],
                "altRefDonorZ": refSpliceInfo["altZScore"],
                "refRefDonorSeq": refSpliceInfo["refSeq"],
                "altRefDonorSeq": refSpliceInfo["altSeq"],
                "refDonorVarStart": refSpliceInfo["varStart"],
                "refDonorVarLength": refSpliceInfo["varLength"],
                "refDonorExonStart": refSpliceInfo["exonStart"],
                "refDonorIntronStart": refSpliceInfo["intronStart"],
                "refDeNovoDonorMES": deNovoSpliceInfo["refMaxEntScanScore"],
                "refDeNovoDonorZ": deNovoSpliceInfo["refZScore"],
                "altDeNovoDonorMES": deNovoSpliceInfo["altMaxEntScanScore"],
                "altDeNovoDonorZ": deNovoSpliceInfo["altZScore"],
                "refDeNovoDonorSeq": deNovoSpliceInfo["refSeq"],
                "altDeNovoDonorSeq": deNovoSpliceInfo["altSeq"],
                "deNovoDonorVarStart": deNovoSpliceInfo["varStart"],
                "deNovoDonorVarLength": deNovoSpliceInfo["varLength"],
                "deNovoDonorExonStart": deNovoSpliceInfo["exonStart"],
                "deNovoDonorIntronStart": deNovoSpliceInfo["intronStart"],
                "deNovoDonorGenomicSplicePos": deNovoSpliceInfo["genomicSplicePos"],
                "deNovoDonorTranscriptSplicePos": deNovoSpliceInfo["transcriptSplicePos"],
                "closestDonorGenomicSplicePos": deNovoSpliceInfo["closestGenomicSplicePos"],
                "closestDonorTranscriptSplicePos": deNovoSpliceInfo["closestTranscriptSplicePos"],
                "closestDonorRefMES": deNovoSpliceInfo["closestRefMaxEntScanScore"],
                "closestDonorRefZ": deNovoSpliceInfo["closestRefZScore"],
                "closestDonorRefSeq": deNovoSpliceInfo["closestRefSeq"],
                "closestDonorAltMES": deNovoSpliceInfo["closestAltMaxEntScanScore"],
                "closestDonorAltZ": deNovoSpliceInfo["closestAltZScore"],
                "closestDonorAltSeq": deNovoSpliceInfo["closestAltSeq"],
                "closestDonorExonStart": deNovoSpliceInfo["closestExonStart"],
                "closestDonorIntronStart": deNovoSpliceInfo["closestIntronStart"],
                "deNovoDonorAltGreaterRefFlag": deNovoSpliceInfo["altGreaterRefFlag"],
                "deNovoDonorAltGreaterClosestRefFlag": deNovoSpliceInfo["altGreaterClosestRefFlag"],
                "deNovoDonorAltGreaterClosestAltFlag": deNovoSpliceInfo["altGreaterClosestAltFlag"],
                "deNovoDonorFrameshiftFlag": deNovoSpliceInfo["frameshiftFlag"],
                "refAccPrior": "N/A",
                "deNovoAccPrior": "N/A",
                "refRefAccMES": "N/A",
                "refRefAccZ": "N/A",
                "altRefAccMES": "N/A",
                "altRefAccZ": "N/A",
                "refRefAccSeq": "N/A",
                "altRefAccSeq": "N/A",
                "refAccVarStart": "N/A",
                "refAccVarLength": "N/A",
                "refAccExonStart": "N/A",
                "refAccIntronStart": "N/A",
                "refDeNovoAccMES": "N/A",
                "refDeNovoAccZ": "N/A",
                "altDeNovoAccMES": "N/A",
                "altDeNovoAccZ": "N/A",
                "refDeNovoAccSeq": "N/A",
                "altDeNovoAccSeq": "N/A",
                "deNovoAccVarStart": "N/A",
                "deNovoAccVarLength": "N/A",
                "deNovoAccExonStart": "N/A",
                "deNovoAccIntronStart": "N/A",
                "deNovoAccGenomicSplicePos": "N/A",
                "deNovoAccTranscriptSplicePos": "N/A",
                "closestAccGenomicSplicePos": "N/A",
                "closestAccTranscriptSplicePos": "N/A",
                "closestAccRefMES": "N/A",
                "closestAccRefZ": "N/A",
                "closestAccRefSeq": "N/A",
                "closestAccAltMES": "N/A",
                "closestAccAltZ": "N/A",
                "closestAccAltSeq": "N/A",
                "closestAccExonStart": "N/A",
                "closestAccIntronStart": "N/A",
                "deNovoAccAltGreaterRefFlag": "N/A",
                "deNovoAccAltGreaterClosestRefFlag": "N/A",
                "deNovoAccAltGreaterClosestAltFlag": "N/A",
                "deNovoAccFrameshiftFlag": "N/A",
                "spliceSite": refSpliceInfo["spliceSite"],
                "spliceRescue": spliceRescue,
                "spliceFlag": spliceFlag,
                "frameshiftFlag": frameshiftFlag,
                "inExonicPortionFlag": inExonicPortionFlag,
                "CIDomainInRegionFlag": CIDomainInRegionFlag,
                "isDivisibleFlag": isDivisibleFlag,
                "lowMESFlag": lowMESFlag,
                "varConsequences": ",".join(varCons)}


def getPriorProbSpliceAcceptorSNS(variant, boundaries, variantData, genome, transcript):
    """
    Given a variant, boundaries (either PRIORS or ENIGMA), and list of dictionaries with variant data
    Determines reference and de novo acceptor scores for variant
      If variant in exon, also determines de novo donor scores and protein prior
    If variant causes a nonsense mutation, determines if splice rescue occurs
    Returns dicitionary containing scores for and sequences ref and de novo splice donor/acceptor,
        closest splice site scores and sequences and protein prior if variant in exon
        "N/A" if score, sequence, or flag not applicable for variant
    Also contains other values:
        applicable prior, highest prior if variant has multiple priors
        applicable classe, highest predicted qualitative enigma class based on applicable prior
        ref prior, prior prob for reference splice sequence
        de novo donor and acceptor priors, prior prob for de novo splice sequence
        splice flag = 1, because variant in reference splice site
        deNovoAccFrameshiftFlag = 1 if de novo acceptor causes a frameshift
        deNovoAccAltGreaterRefFlag = 1 if variant is possible de novo acceptor variant
        deNovoAccAltGreaterClosestRefFlag = 1 if variant is possible de novo acceptor variant
        deNovoAccAltGreaterClosestAltFlag = 1 if variant is possible de novo acceptor variant
        deNovoDonroFrameshiftFlag =1 if de novo donor causes a frameshift
        deNovoDonorAltGreaterRefFlag = 1 if variant is possible de novo donor variant
        deNovoDonorAltGreaterClosestRefFlag = 1 if variant is possible de novo donor variant
        deNovoDonorAltGreaterClosestAltFlag = 1 if variant is possible de novo donor variant
        spliceRescue = 1 if splice rescue possible for nonsense variant, 0 otherwise, N/A if not nonsense variant
        spliceFlag = 1 if splice rescue is possible so variant can be flagged for further analysis, 0 otherwise
        frameshiftFlag = 1 if nonsense variant causes a frameshift mutation also, 0 if not, N/A if not nonsense variant
        inExonicPortionFlag = 1 if variant IS in exonic portion of highest scoring window and variant is a nonsense variant,
           0 if NOT in exonic portion, N/A if not nonsense variant
        CIDomainRegionFlag = 1 if truncating region (between new donor and next reference splice acceptor)
           includes a clinically import domain for nonsense variant, 0 if not, N/A if not nonsense variant
        isDivisibleFlag = 1 if distance between de novo donor and wild-type donor splice site is NOT divisible by 3 for nonsense variant,
           0 if it IS divisible, N/A if not a nonsense variant
    Dictionary also contains formatting variables for each listed sequence
    """
    if compute.varInSpliceRegion(variant, donor=False, deNovo=False) and extract.getVarType(variant) == "substitution":
        refSpliceInfo = getPriorProbRefSpliceAcceptorSNS(variant, boundaries)
        deNovoAccInfo = getPriorProbDeNovoAcceptorSNS(variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH, genome,
                                                      transcript)
        refPrior = refSpliceInfo["priorProb"]
        proteinPrior = "N/A"
        applicablePrior = refSpliceInfo["priorProb"]
        if compute.varInExon(variant):
            deNovoDonorInfo = getPriorProbDeNovoDonorSNS(variant, boundaries, STD_EXONIC_PORTION, genome, transcript,
                                                         deNovoDonorInRefAcc=True)
            deNovoDonorPrior = deNovoDonorInfo["priorProb"]
            proteinInfo = getPriorProbProteinSNS(variant, variantData)
            proteinPrior = proteinInfo["priorProb"]
            if deNovoDonorPrior != "N/A" and proteinPrior != "N/A":
                applicablePrior = max(deNovoDonorPrior, proteinPrior, refPrior)
            else:
                applicablePrior = max(proteinPrior, refPrior)
        else:
            deNovoDonorPrior = "N/A"
            deNovoDonorInfo = {"refMaxEntScanScore": "N/A",
                               "refZScore": "N/A",
                               "altMaxEntScanScore": "N/A",
                               "altZScore": "N/A",
                               "refSeq": "N/A",
                               "altSeq": "N/A",
                               "varStart": "N/A",
                               "varLength": "N/A",
                               "exonStart": "N/A",
                               "intronStart": "N/A",
                               "genomicSplicePos": "N/A",
                               "transcriptSplicePos": "N/A",
                               "closestGenomicSplicePos": "N/A",
                               "closestTranscriptSplicePos": "N/A",
                               "closestRefMaxEntScanScore": "N/A",
                               "closestRefZScore": "N/A",
                               "closestRefSeq": "N/A",
                               "closestAltMaxEntScanScore": "N/A",
                               "closestAltZScore": "N/A",
                               "closestAltSeq": "N/A",
                               "closestExonStart": "N/A",
                               "closestIntronStart": "N/A",
                               "altGreaterRefFlag": "N/A",
                               "altGreaterClosestRefFlag": "N/A",
                               "altGreaterClosestAltFlag": "N/A",
                               "frameshiftFlag": "N/A",
                               "priorProb": "N/A",
                               "enigmaClass": "N/A"}

        spliceRescue = "N/A"
        spliceFlag = 0
        frameshiftFlag = "N/A"
        inExonicPortionFlag = "N/A"
        CIDomainInRegionFlag = "N/A"
        isDivisibleFlag = "N/A"
        lowMESFlag = "N/A"
        # to check for nonsense variants in exonic portion of splice acceptor site
        varCons = extract.getVarConsequences(variant)
        if compute.varInExon(variant) and "stop_gained" in varCons:
            nonsenseData = getPriorProbSpliceRescueNonsenseSNS(variant, boundaries, deNovoDonorInRefAcc=True)
            applicablePrior = nonsenseData["priorProb"]
            spliceRescue = nonsenseData["spliceRescue"]
            spliceFlag = nonsenseData["spliceFlag"]
            frameshiftFlag = nonsenseData["frameshiftFlag"]
            inExonicPortionFlag = nonsenseData["inExonicPortionFlag"]
            CIDomainInRegionFlag = nonsenseData["CIDomainInRegionFlag"]
            isDivisibleFlag = nonsenseData["isDivisibleFlag"]
            lowMESFlag = nonsenseData["lowMESFlag"]

        return {"applicablePrior": applicablePrior,
                "applicableEnigmaClass": extract.getEnigmaClass(applicablePrior),
                "proteinPrior": proteinPrior,
                "refDonorPrior": "N/A",
                "deNovoDonorPrior": deNovoDonorPrior,
                "refRefDonorMES": "N/A",
                "refRefDonorZ": "N/A",
                "altRefDonorMES": "N/A",
                "altRefDonorZ": "N/A",
                "refRefDonorSeq": "N/A",
                "altRefDonorSeq": "N/A",
                "refDonorVarStart": "N/A",
                "refDonorVarLength": "N/A",
                "refDonorExonStart": "N/A",
                "refDonorIntronStart": "N/A",
                "refDeNovoDonorMES": deNovoDonorInfo["refMaxEntScanScore"],
                "refDeNovoDonorZ": deNovoDonorInfo["refZScore"],
                "altDeNovoDonorMES": deNovoDonorInfo["altMaxEntScanScore"],
                "altDeNovoDonorZ": deNovoDonorInfo["altZScore"],
                "refDeNovoDonorSeq": deNovoDonorInfo["refSeq"],
                "altDeNovoDonorSeq": deNovoDonorInfo["altSeq"],
                "deNovoDonorVarStart": deNovoDonorInfo["varStart"],
                "deNovoDonorVarLength": deNovoDonorInfo["varLength"],
                "deNovoDonorExonStart": deNovoDonorInfo["exonStart"],
                "deNovoDonorIntronStart": deNovoDonorInfo["intronStart"],
                "deNovoDonorGenomicSplicePos": deNovoDonorInfo["genomicSplicePos"],
                "deNovoDonorTranscriptSplicePos": deNovoDonorInfo["transcriptSplicePos"],
                "closestDonorGenomicSplicePos": deNovoDonorInfo["closestGenomicSplicePos"],
                "closestDonorTranscriptSplicePos": deNovoDonorInfo["closestTranscriptSplicePos"],
                "closestDonorRefMES": deNovoDonorInfo["closestRefMaxEntScanScore"],
                "closestDonorRefZ": deNovoDonorInfo["closestRefZScore"],
                "closestDonorRefSeq": deNovoDonorInfo["closestRefSeq"],
                "closestDonorAltMES": deNovoDonorInfo["closestAltMaxEntScanScore"],
                "closestDonorAltZ": deNovoDonorInfo["closestAltZScore"],
                "closestDonorAltSeq": deNovoDonorInfo["closestAltSeq"],
                "closestDonorExonStart": deNovoDonorInfo["closestExonStart"],
                "closestDonorIntronStart": deNovoDonorInfo["closestIntronStart"],
                "deNovoDonorAltGreaterRefFlag": deNovoDonorInfo["altGreaterRefFlag"],
                "deNovoDonorAltGreaterClosestRefFlag": deNovoDonorInfo["altGreaterClosestRefFlag"],
                "deNovoDonorAltGreaterClosestAltFlag": deNovoDonorInfo["altGreaterClosestAltFlag"],
                "deNovoDonorFrameshiftFlag": deNovoDonorInfo["frameshiftFlag"],
                "deNovoAccPrior": deNovoAccInfo["priorProb"],
                "refAccPrior": refPrior,
                "refRefAccMES": refSpliceInfo["refMaxEntScanScore"],
                "refRefAccZ": refSpliceInfo["refZScore"],
                "altRefAccMES": refSpliceInfo["altMaxEntScanScore"],
                "altRefAccZ": refSpliceInfo["altZScore"],
                "refRefAccSeq": refSpliceInfo["refSeq"],
                "altRefAccSeq": refSpliceInfo["altSeq"],
                "refAccVarStart": refSpliceInfo["varStart"],
                "refAccVarLength": refSpliceInfo["varLength"],
                "refAccExonStart": refSpliceInfo["exonStart"],
                "refAccIntronStart": refSpliceInfo["intronStart"],
                "refDeNovoAccMES": deNovoAccInfo["refMaxEntScanScore"],
                "refDeNovoAccZ": deNovoAccInfo["refZScore"],
                "altDeNovoAccMES": deNovoAccInfo["altMaxEntScanScore"],
                "altDeNovoAccZ": deNovoAccInfo["altZScore"],
                "refDeNovoAccSeq": deNovoAccInfo["refSeq"],
                "altDeNovoAccSeq": deNovoAccInfo["altSeq"],
                "deNovoAccVarStart": deNovoAccInfo["varStart"],
                "deNovoAccVarLength": deNovoAccInfo["varLength"],
                "deNovoAccExonStart": deNovoAccInfo["exonStart"],
                "deNovoAccIntronStart": deNovoAccInfo["intronStart"],
                "deNovoAccGenomicSplicePos": deNovoAccInfo["genomicSplicePos"],
                "deNovoAccTranscriptSplicePos": deNovoAccInfo["transcriptSplicePos"],
                "closestAccGenomicSplicePos": deNovoAccInfo["closestGenomicSplicePos"],
                "closestAccTranscriptSplicePos": deNovoAccInfo["closestTranscriptSplicePos"],
                "closestAccRefMES": deNovoAccInfo["closestRefMaxEntScanScore"],
                "closestAccRefZ": deNovoAccInfo["closestRefZScore"],
                "closestAccRefSeq": deNovoAccInfo["closestRefSeq"],
                "closestAccAltMES": deNovoAccInfo["closestAltMaxEntScanScore"],
                "closestAccAltZ": deNovoAccInfo["closestAltZScore"],
                "closestAccAltSeq": deNovoAccInfo["closestAltSeq"],
                "closestAccExonStart": deNovoAccInfo["closestExonStart"],
                "closestAccIntronStart": deNovoAccInfo["closestIntronStart"],
                "deNovoAccAltGreaterRefFlag": deNovoAccInfo["altGreaterRefFlag"],
                "deNovoAccAltGreaterClosestRefFlag": deNovoAccInfo["altGreaterClosestRefFlag"],
                "deNovoAccAltGreaterClosestAltFlag": deNovoAccInfo["altGreaterClosestAltFlag"],
                "deNovoAccFrameshiftFlag": deNovoAccInfo["frameshiftFlag"],
                "spliceSite": refSpliceInfo["spliceSite"],
                "spliceRescue": spliceRescue,
                "spliceFlag": spliceFlag,
                "frameshiftFlag": frameshiftFlag,
                "inExonicPortionFlag": inExonicPortionFlag,
                "CIDomainInRegionFlag": CIDomainInRegionFlag,
                "isDivisibleFlag": isDivisibleFlag,
                "lowMESFlag": lowMESFlag,
                "varConsequences": ",".join(varCons)}


def getPriorProbProteinSNS(variant, variantData):
    """
    Given a variant and a list of dictionaries containing variant data,
    Returns a dictionary containing:
      the variant's protein prior probability and enigma class for that prior
    """
    proteinPrior = "-"
    enigmaClass = "-"
    if extract.getVarType(variant) == "substitution":
        varHGVS = variant["HGVS_cDNA"]
        if varHGVS == "-":
            # if HGVS_cDNA field is blank use pyhgvs field instead
            # [12:] parses out the NM_ accession so varHGVS is in format c.65C>T
            varHGVS = variant["pyhgvs_cDNA"][12:]
        varGene = variant["Gene_Symbol"]

        foundVar = variantData[(varGene, varHGVS)]
        proteinPrior = float(foundVar["protein_prior"])
        enigmaClass = extract.getEnigmaClass(proteinPrior)

        return {"priorProb": proteinPrior,
                "enigmaClass": enigmaClass}


def getPriorProbInGreyZoneSNS(variant, boundaries, variantData):
    """
    Given a variant and a list of dicitionaries with variant data,
    Returns applicable prior and enigma class based on protein priors for that variant
    Dictionary also contains other values that are either "N/A", "-", or 0 because they are not relevant
    """
    if extract.getVarType(variant) == "substitution" and compute.getVarLocation(variant, boundaries) == "grey_zone_variant":
        proteinData = getPriorProbProteinSNS(variant, variantData)
        proteinPrior = proteinData["priorProb"]
        if proteinPrior == PATHOGENIC_PROBABILITY:
            proteinPrior = CAPPED_PROBABILITY

        return {"applicablePrior": proteinPrior,
                "applicableEnigmaClass": extract.getEnigmaClass(proteinPrior),
                "proteinPrior": proteinPrior,
                "refDonorPrior": "N/A",
                "deNovoDonorPrior": "N/A",
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
                "deNovoDonorAltGreaterRefFlag": "N/A",
                "deNovoDonorAltGreaterClosestRefFlag": "N/A",
                "deNovoDonorAltGreaterClosestAltFlag": "N/A",
                "deNovoDonorFrameshiftFlag": "N/A",
                "refAccPrior": "N/A",
                "deNovoAccPrior": "N/A",
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
                "deNovoAccAltGreaterRefFlag": "N/A",
                "deNovoAccAltGreaterClosestRefFlag": "N/A",
                "deNovoAccAltGreaterClosestAltFlag": "N/A",
                "deNovoAccFrameshiftFlag": "N/A",
                "spliceSite": 0,
                "spliceRescue": "N/A",
                "spliceFlag": 0,
                "frameshiftFlag": "N/A",
                "inExonicPortionFlag": "N/A",
                "CIDomainInRegionFlag": "N/A",
                "isDivisibleFlag": "N/A",
                "lowMESFlag": "N/A"}


def getPriorProbInExonSNS(variant, boundaries, variantData, genome, transcript):
    """
    Given a variant, boundaries (either "enigma" or "priors") and a list of dictionaries containing variant data:
      1. Checks that variant is in an exon or clinically important domains and NOT in a splice site
      2. Checks that variant is a SNS variant
      3. Gets protein prior from variantData
      4. Determines if variant is a nonsense variant, if yes determines if splice rescue occurs
      5. If not a nonsense variant, calculates de novo donor prior and de novo acceptor prior if applicable
         Gets applicable prior if variant has a de novo donor prior
    Genome is a SequenceFileDB for genome and transcript is a pyhgvs transcript object)
       both genome and transcript are necessary to convert from genomic to transcript coordinates
    Returns a dictionary containing all values, dictionary entry is "-" if not relevant to variant
    Values in dictionary include:
        applicable prior, highest prior if variant has multiple priors
        applicable class, highest predicted qualitative enigma class
        splice site = 0 because these variants are not in reference splice sites
        de novo donor and acceptor priors, prior prob for de novo splice sequence
        de novo acc flags (altGreaterRef, altGreaterClosestRef, altGreaterClosestAlt) = 1 if variant is possible de novo acceptor variant
        de novo donor flags (altGreaterRef, altGreaterClosestRef, altGreaterClosestAlt) = 1 if variant is possible de novo donor variant
        de novo acc/donor frameshift flags = 1 if de novo donor or acceptor causes a frameshfit
        spliceRescue = 1 if splice rescue possible for nonsense variant, 0 otherwise
        spliceFlag = 1 if splice rescue is possible so variant can be flagged for further analysis, 0 otherwise
        frameshiftFlag = 1 if nonsense variant causes a frameshift mutation also, 0 if not, N/A if not nonsense variant
        inExonicPortionFlag = 1 if variant IS in exonic portion of highest scoring window and variant is a nonsense variant,
           0 if NOT in exonic portion, N/A if not nonsense variant
        CIDomainRegionFlag = 1 if truncating region (between new donor and next reference splice acceptor)
           includes a clinically import domain for nonsense variant, 0 if not, N/A if not nonsense variant
        isDivisibleFlag = 1 if distance between de novo donor and wild-type donor splice site is NOT divisible by 3 for nonsense variant,
           0 if it IS divisible, N/A if not a nonsense variant
    """
    varLoc = compute.getVarLocation(variant, boundaries)
    if (varLoc == "exon_variant" or "CI_domain_variant") and extract.getVarType(variant) == "substitution":
        proteinData = getPriorProbProteinSNS(variant, variantData)
        spliceRescue = "N/A"
        spliceFlag = 0
        frameshiftFlag = "N/A"
        inExonicPortionFlag = "N/A"
        CIDomainInRegionFlag = "N/A"
        isDivisibleFlag = "N/A"
        lowMESFlag = "N/A"
        deNovoDonorData = getPriorProbDeNovoDonorSNS(variant, boundaries, STD_EXONIC_PORTION, genome, transcript,
                                                     deNovoDonorInRefAcc=False)
        if compute.varInSpliceRegion(variant, donor=False, deNovo=True):
            deNovoAccData = getPriorProbDeNovoAcceptorSNS(variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH, genome,
                                                          transcript)
        else:
            deNovoAccData = {"priorProb": "N/A",
                             "refMaxEntScanScore": "N/A",
                             "altMaxEntScanScore": "N/A",
                             "refZScore": "N/A",
                             "altZScore": "N/A",
                             "refSeq": "N/A",
                             "altSeq": "N/A",
                             "varStart": "N/A",
                             "varLength": "N/A",
                             "exonStart": "N/A",
                             "intronStart": "N/A",
                             "genomicSplicePos": "N/A",
                             "transcriptSplicePos": "N/A",
                             "closestGenomicSplicePos": "N/A",
                             "closestTranscriptSplicePos": "N/A",
                             "closestRefMaxEntScanScore": "N/A",
                             "closestRefZScore": "N/A",
                             "closestRefSeq": "N/A",
                             "closestAltMaxEntScanScore": "N/A",
                             "closestAltZScore": "N/A",
                             "closestAltSeq": "N/A",
                             "closestExonStart": "N/A",
                             "closestIntronStart": "N/A",
                             "altGreaterRefFlag": "N/A",
                             "altGreaterClosestRefFlag": "N/A",
                             "altGreaterClosestAltFlag": "N/A",
                             "frameshiftFlag": "N/A"}
        varCons = extract.getVarConsequences(variant)
        if "stop_gained" in varCons:
            nonsenseData = getPriorProbSpliceRescueNonsenseSNS(variant, boundaries, deNovoDonorInRefAcc=False)
            applicablePrior = nonsenseData["priorProb"]
            applicableClass = nonsenseData["enigmaClass"]
            spliceRescue = nonsenseData["spliceRescue"]
            spliceFlag = nonsenseData["spliceFlag"]
            frameshiftFlag = nonsenseData["frameshiftFlag"]
            inExonicPortionFlag = nonsenseData["inExonicPortionFlag"]
            CIDomainInRegionFlag = nonsenseData["CIDomainInRegionFlag"]
            isDivisibleFlag = nonsenseData["isDivisibleFlag"]
            lowMESFlag = nonsenseData["lowMESFlag"]
        else:
            applicablePrior = proteinData["priorProb"]
            applicableClass = proteinData["enigmaClass"]
            if deNovoDonorData["priorProb"] != "N/A":
                applicablePrior = max(proteinData["priorProb"], deNovoDonorData["priorProb"])
                applicableClass = extract.getEnigmaClass(applicablePrior)

        return {"applicablePrior": applicablePrior,
                "applicableEnigmaClass": applicableClass,
                "proteinPrior": proteinData["priorProb"],
                "refDonorPrior": "N/A",
                "deNovoDonorPrior": deNovoDonorData["priorProb"],
                "refRefDonorMES": "N/A",
                "refRefDonorZ": "N/A",
                "altRefDonorMES": "N/A",
                "altRefDonorZ": "N/A",
                "refRefDonorSeq": "N/A",
                "altRefDonorSeq": "N/A",
                "refDonorVarStart": "N/A",
                "refDonorVarLength": "N/A",
                "refDonorExonStart": "N/A",
                "refDonorIntronStart": "N/A",
                "refDeNovoDonorMES": deNovoDonorData["refMaxEntScanScore"],
                "refDeNovoDonorZ": deNovoDonorData["refZScore"],
                "altDeNovoDonorMES": deNovoDonorData["altMaxEntScanScore"],
                "altDeNovoDonorZ": deNovoDonorData["altZScore"],
                "refDeNovoDonorSeq": deNovoDonorData["refSeq"],
                "altDeNovoDonorSeq": deNovoDonorData["altSeq"],
                "deNovoDonorVarStart": deNovoDonorData["varStart"],
                "deNovoDonorVarLength": deNovoDonorData["varLength"],
                "deNovoDonorExonStart": deNovoDonorData["exonStart"],
                "deNovoDonorIntronStart": deNovoDonorData["intronStart"],
                "deNovoDonorGenomicSplicePos": deNovoDonorData["genomicSplicePos"],
                "deNovoDonorTranscriptSplicePos": deNovoDonorData["transcriptSplicePos"],
                "closestDonorGenomicSplicePos": deNovoDonorData["closestGenomicSplicePos"],
                "closestDonorTranscriptSplicePos": deNovoDonorData["closestTranscriptSplicePos"],
                "closestDonorRefMES": deNovoDonorData["closestRefMaxEntScanScore"],
                "closestDonorRefZ": deNovoDonorData["closestRefZScore"],
                "closestDonorRefSeq": deNovoDonorData["closestRefSeq"],
                "closestDonorAltMES": deNovoDonorData["closestAltMaxEntScanScore"],
                "closestDonorAltZ": deNovoDonorData["closestAltZScore"],
                "closestDonorAltSeq": deNovoDonorData["closestAltSeq"],
                "closestDonorExonStart": deNovoDonorData["closestExonStart"],
                "closestDonorIntronStart": deNovoDonorData["closestIntronStart"],
                "deNovoDonorAltGreaterRefFlag": deNovoDonorData["altGreaterRefFlag"],
                "deNovoDonorAltGreaterClosestRefFlag": deNovoDonorData["altGreaterClosestRefFlag"],
                "deNovoDonorAltGreaterClosestAltFlag": deNovoDonorData["altGreaterClosestAltFlag"],
                "deNovoDonorFrameshiftFlag": deNovoDonorData["frameshiftFlag"],
                "deNovoAccPrior": deNovoAccData["priorProb"],
                "refAccPrior": "N/A",
                "refRefAccMES": "N/A",
                "refRefAccZ": "N/A",
                "altRefAccMES": "N/A",
                "altRefAccZ": "N/A",
                "refRefAccSeq": "N/A",
                "altRefAccSeq": "N/A",
                "refAccVarStart": "N/A",
                "refAccVarLength": "N/A",
                "refAccExonStart": "N/A",
                "refAccIntronStart": "N/A",
                "refDeNovoAccMES": deNovoAccData["refMaxEntScanScore"],
                "refDeNovoAccZ": deNovoAccData["refZScore"],
                "altDeNovoAccMES": deNovoAccData["altMaxEntScanScore"],
                "altDeNovoAccZ": deNovoAccData["altZScore"],
                "refDeNovoAccSeq": deNovoAccData["refSeq"],
                "altDeNovoAccSeq": deNovoAccData["altSeq"],
                "deNovoAccVarStart": deNovoAccData["varStart"],
                "deNovoAccVarLength": deNovoAccData["varLength"],
                "deNovoAccExonStart": deNovoAccData["exonStart"],
                "deNovoAccIntronStart": deNovoAccData["intronStart"],
                "deNovoAccGenomicSplicePos": deNovoAccData["genomicSplicePos"],
                "deNovoAccTranscriptSplicePos": deNovoAccData["transcriptSplicePos"],
                "closestAccGenomicSplicePos": deNovoAccData["closestGenomicSplicePos"],
                "closestAccTranscriptSplicePos": deNovoAccData["closestTranscriptSplicePos"],
                "closestAccRefMES": deNovoAccData["closestRefMaxEntScanScore"],
                "closestAccRefZ": deNovoAccData["closestRefZScore"],
                "closestAccRefSeq": deNovoAccData["closestRefSeq"],
                "closestAccAltMES": deNovoAccData["closestAltMaxEntScanScore"],
                "closestAccAltZ": deNovoAccData["closestAltZScore"],
                "closestAccAltSeq": deNovoAccData["closestAltSeq"],
                "closestAccExonStart": deNovoAccData["closestExonStart"],
                "closestAccIntronStart": deNovoAccData["closestIntronStart"],
                "deNovoAccAltGreaterRefFlag": deNovoAccData["altGreaterRefFlag"],
                "deNovoAccAltGreaterClosestRefFlag": deNovoAccData["altGreaterClosestRefFlag"],
                "deNovoAccAltGreaterClosestAltFlag": deNovoAccData["altGreaterClosestAltFlag"],
                "deNovoAccFrameshiftFlag": deNovoAccData["frameshiftFlag"],
                "spliceSite": 0,
                "spliceRescue": spliceRescue,
                "spliceFlag": spliceFlag,
                "frameshiftFlag": frameshiftFlag,
                "inExonicPortionFlag": inExonicPortionFlag,
                "CIDomainInRegionFlag": CIDomainInRegionFlag,
                "isDivisibleFlag": isDivisibleFlag,
                "lowMESFlag": lowMESFlag,
                "varConsequences": ",".join(varCons)}


def getPriorProbOutsideTranscriptBoundsSNS(variant, boundaries):
    """
    Given a variant and boundaries (either "enigma" or "priors"),
    Checks that variant is outside transcript boundaries
    Returns prior prob and predicted qualitative enigma class
    Dictionary also contains other values that are either "N/A", "-", or 0 because they are not relevant
    """
    varLoc = compute.getVarLocation(variant, boundaries)
    varType = extract.getVarType(variant)
    if varLoc == "outside_transcript_boundaries_variant" and varType == "substitution":
        priorProb = LOW_PROBABILITY
        return {"applicablePrior": priorProb,
                "applicableEnigmaClass": extract.getEnigmaClass(priorProb),
                "proteinPrior": "N/A",
                "refDonorPrior": "N/A",
                "deNovoDonorPrior": "N/A",
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
                "deNovoDonorAltGreaterRefFlag": "N/A",
                "deNovoDonorAltGreaterClosestRefFlag": "N/A",
                "deNovoDonorAltGreaterClosestAltFlag": "N/A",
                "deNovoDonorFrameshiftFlag": "N/A",
                "refAccPrior": "N/A",
                "deNovoAccPrior": "N/A",
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
                "deNovoAccAltGreaterRefFlag": "N/A",
                "deNovoAccAltGreaterClosestRefFlag": "N/A",
                "deNovoAccAltGreaterClosestAltFlag": "N/A",
                "deNovoAccFrameshiftFlag": "N/A",
                "spliceSite": 0,
                "spliceRescue": "N/A",
                "spliceFlag": 0,
                "frameshiftFlag": "N/A",
                "inExonicPortionFlag": "N/A",
                "CIDomainInRegionFlag": "N/A",
                "isDivisibleFlag": "N/A",
                "lowMESFlag": "N/A"}


def getPriorProbIntronicDeNovoDonorSNS(variant, genome, transcript):
    """
    Given a variant,
      1. Checks that variant is NOT in exon or reference donor/acceptor site
      2. Checks that variant is a substitution variant
    Genome is a SequenceFileDB for genome and transcript is a pyhgvs transcript object)
       both genome and transcript are necessary to convert from genomic to transcript coordinates
    Determines if alt MES score is greater than ref MES score for highest scoring sliding window
      If altMES > refMES, then deNovoDonorAltGreaterRefFlag = 1 (0 otherwise)
    Determines closest ref donor scores
      If altMES > closestRefDonorMES, then deNovoDonorAltGreaterClosestRefFlag = 1 (0 otherwise)
    If either flag is equal to 1, flag variant for further analysis (spliceFlag = 1), spliceFlag = 0 otherwise
    altGreaterClosestAltFlag is always equals N/A because this function is not used for any variants in ref splice sites
    Returns dictionary containing prior prob, enigma class, de novo donor scores, and splice flag
    Also contains closest ref donor scores and de novo and closest donor genomic and transcript splice positions
    """
    inExon = compute.varInExon(variant)
    inRefDonor = compute.varInSpliceRegion(variant, donor=True, deNovo=False)
    inRefAcc = compute.varInSpliceRegion(variant, donor=False, deNovo=False)
    if not inExon and not inRefDonor and not inRefAcc:
        if extract.getVarType(variant) == "substitution":
            deNovoDonorInfo = compute.getMaxMaxEntScanScoreSlidingWindowSNS(variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH,
                                                                    donor=True, deNovo=False, deNovoDonorInRefAcc=False)
            # +-1 because per HCI PRIORS website donor position is defined as being first nucleotide that is NOT included in spliced exon
            if verify.getVarStrand(variant) == "-":
                newGenomicSplicePos = extract.getNewSplicePosition(variant["Pos"], "-", deNovoDonorInfo["varWindowPosition"],
                                                           deNovoDonorInfo["inExonicPortion"], STD_EXONIC_PORTION,
                                                           STD_ACC_INTRONIC_LENGTH,
                                                           donor=True) - 1
            else:
                newGenomicSplicePos = extract.getNewSplicePosition(variant["Pos"], "+", deNovoDonorInfo["varWindowPosition"],
                                                           deNovoDonorInfo["inExonicPortion"], STD_EXONIC_PORTION,
                                                           STD_ACC_INTRONIC_LENGTH,
                                                           donor=True) + 1
            refMES = deNovoDonorInfo["refMaxEntScanScore"]
            altMES = deNovoDonorInfo["altMaxEntScanScore"]
            deNovoOffset = 0
            closestDonorInfo = compute.getClosestSpliceSiteScores(variant, deNovoOffset, donor=True, deNovo=False,
                                                          deNovoDonorInRefAcc=False, testMode=False)
            closestMES = closestDonorInfo["maxEntScanScore"]
            spliceFlag = 0
            altGreaterRefFlag = 0
            altGreaterClosestRefFlag = 0
            if altMES > refMES:
                spliceFlag = 1
                altGreaterRefFlag = 1
            if altMES > closestMES:
                spliceFlag = 1
                altGreaterClosestRefFlag = 1

            frameshiftStatus = compute.getDeNovoSpliceFrameshiftStatus(variant, donor=True, deNovoDonorInRefAcc=False)
            frameshiftFlag = 0
            if frameshiftStatus:
                frameshiftFlag = 1

            # converts genomic splice position to transcript splice position
            newTranscriptSplicePos = verify.convertGenomicPosToTranscriptPos(newGenomicSplicePos, extract.getVarChrom(variant), genome,
                                                                      transcript)
            # converts closest genomic splice position to transcript splice position
            closestGenomicSplicePos = closestDonorInfo["genomicSplicePos"]
            closestTranscriptSplicePos = verify.convertGenomicPosToTranscriptPos(closestGenomicSplicePos, extract.getVarChrom(variant),
                                                                          genome, transcript)

            return {"priorProb": "N/A",
                    "enigmaClass": "N/A",
                    "refMaxEntScanScore": refMES,
                    "refZScore": deNovoDonorInfo["refZScore"],
                    "altMaxEntScanScore": altMES,
                    "altZScore": deNovoDonorInfo["altZScore"],
                    "refSeq": deNovoDonorInfo["refSeq"],
                    "altSeq": deNovoDonorInfo["altSeq"],
                    "varStart": deNovoDonorInfo["varStart"],
                    "varLength": deNovoDonorInfo["varLength"],
                    "exonStart": 0,
                    "intronStart": STD_EXONIC_PORTION,
                    "genomicSplicePos": verify.formatSplicePosition(newGenomicSplicePos, transcript=False),
                    "transcriptSplicePos": verify.formatSplicePosition(newTranscriptSplicePos, transcript=True),
                    "closestGenomicSplicePos": verify.formatSplicePosition(closestGenomicSplicePos, transcript=False),
                    "closestTranscriptSplicePos": verify.formatSplicePosition(closestTranscriptSplicePos, transcript=True),
                    "closestRefSeq": closestDonorInfo["sequence"],
                    "closestRefMaxEntScanScore": closestMES,
                    "closestRefZScore": closestDonorInfo["zScore"],
                    "closestAltSeq": "N/A",
                    "closestAltMaxEntScanScore": "N/A",
                    "closestAltZScore": "N/A",
                    "closestExonStart": closestDonorInfo["exonStart"],
                    "closestIntronStart": closestDonorInfo["intronStart"],
                    "altGreaterRefFlag": altGreaterRefFlag,
                    "altGreaterClosestRefFlag": altGreaterClosestRefFlag,
                    "altGreaterClosestAltFlag": "N/A",
                    "frameshiftFlag": frameshiftFlag,
                    "spliceFlag": spliceFlag}


def getPriorProbInIntronSNS(variant, boundaries, genome, transcript):
    """
    Given a variant and boundaries (either "priors or "enigma"),
    Checks that variant is located in an intron and is a substitution variant
    Determines if variant creates a de novo donor site in the intron
    Genome is a SequenceFileDB for genome and transcript is a pyhgvs transcript object)
       both genome and transcript are necessary to convert from genomic to transcript coordinates
    Returns a dictionary containing applicable prior and predicted qualitative enigma class
    Dictionary also contains de novo donor ref and alt scores
    AND a spliceFlag which is equal to 1 if variant creates a better de novo splice site than ref sequence
    Rest of values in dictionary are equal to 0, "-", or N/A because they are not relevant to variant
    """
    varLoc = compute.getVarLocation(variant, boundaries)
    varType = extract.getVarType(variant)
    if varLoc == "intron_variant" and varType == "substitution":
        deNovoDonorData = getPriorProbIntronicDeNovoDonorSNS(variant, genome, transcript)
        if deNovoDonorData["spliceFlag"] == 0:
            priorProb = LOW_PROBABILITY
            enigmaClass = extract.getEnigmaClass(priorProb)
        else:
            priorProb = deNovoDonorData["priorProb"]
            enigmaClass = deNovoDonorData["enigmaClass"]

        return {"applicablePrior": priorProb,
                "applicableEnigmaClass": enigmaClass,
                "proteinPrior": "N/A",
                "refDonorPrior": "N/A",
                "deNovoDonorPrior": "N/A",
                "refRefDonorMES": "N/A",
                "refRefDonorZ": "N/A",
                "altRefDonorMES": "N/A",
                "altRefDonorZ": "N/A",
                "refRefDonorSeq": "N/A",
                "altRefDonorSeq": "N/A",
                "refDonorVarStart": "N/A",
                "refDonorVarLength": "N/A",
                "refDonorExonStart": "N/A",
                "refDonorIntronStart": "N/A",
                "refDeNovoDonorMES": deNovoDonorData["refMaxEntScanScore"],
                "refDeNovoDonorZ": deNovoDonorData["refZScore"],
                "altDeNovoDonorMES": deNovoDonorData["altMaxEntScanScore"],
                "altDeNovoDonorZ": deNovoDonorData["altZScore"],
                "refDeNovoDonorSeq": deNovoDonorData["refSeq"],
                "altDeNovoDonorSeq": deNovoDonorData["altSeq"],
                "deNovoDonorVarStart": deNovoDonorData["varStart"],
                "deNovoDonorVarLength": deNovoDonorData["varLength"],
                "deNovoDonorExonStart": deNovoDonorData["exonStart"],
                "deNovoDonorIntronStart": deNovoDonorData["intronStart"],
                "deNovoDonorGenomicSplicePos": deNovoDonorData["genomicSplicePos"],
                "deNovoDonorTranscriptSplicePos": deNovoDonorData["transcriptSplicePos"],
                "closestDonorGenomicSplicePos": deNovoDonorData["closestGenomicSplicePos"],
                "closestDonorTranscriptSplicePos": deNovoDonorData["closestTranscriptSplicePos"],
                "closestDonorRefMES": deNovoDonorData["closestRefMaxEntScanScore"],
                "closestDonorRefZ": deNovoDonorData["closestRefZScore"],
                "closestDonorRefSeq": deNovoDonorData["closestRefSeq"],
                "closestDonorAltMES": "N/A",
                "closestDonorAltZ": "N/A",
                "closestDonorAltSeq": "N/A",
                "closestDonorExonStart": deNovoDonorData["closestExonStart"],
                "closestDonorIntronStart": deNovoDonorData["closestIntronStart"],
                "deNovoDonorAltGreaterRefFlag": deNovoDonorData["altGreaterRefFlag"],
                "deNovoDonorAltGreaterClosestRefFlag": deNovoDonorData["altGreaterClosestRefFlag"],
                "deNovoDonorAltGreaterClosestAltFlag": deNovoDonorData["altGreaterClosestAltFlag"],
                "deNovoDonorFrameshiftFlag": deNovoDonorData["frameshiftFlag"],
                "deNovoAccPrior": "N/A",
                "refAccPrior": "N/A",
                "refRefAccMES": "N/A",
                "refRefAccZ": "N/A",
                "altRefAccMES": "N/A",
                "altRefAccZ": "N/A",
                "refRefAccSeq": "N/A",
                "altRefAccSeq": "N/A",
                "refAccVarStart": "N/A",
                "refAccVarLength": "N/A",
                "refAccExonStart": "N/A",
                "refAccIntronStart": "N/A",
                "refDeNovoAccMES": "N/A",
                "refDeNovoAccZ": "N/A",
                "altDeNovoAccMES": "N/A",
                "altDeNovoAccZ": "N/A",
                "refDeNovoAccSeq": "N/A",
                "altDeNovoAccSeq": "N/A",
                "deNovoAccVarStart": "N/A",
                "deNovoAccVarLength": "N/A",
                "deNovoAccExonStart": "N/A",
                "deNovoAccIntronStart": "N/A",
                "deNovoAccGenomicSplicePos": "N/A",
                "deNovoAccTranscriptSplicePos": "N/A",
                "closestAccGenomicSplicePos": "N/A",
                "closestAccTranscriptSplicePos": "N/A",
                "closestAccRefMES": "N/A",
                "closestAccRefZ": "N/A",
                "closestAccRefSeq": "N/A",
                "closestAccAltMES": "N/A",
                "closestAccAltZ": "N/A",
                "closestAccAltSeq": "N/A",
                "closestAccExonStart": "N/A",
                "closestAccIntronStart": "N/A",
                "deNovoAccAltGreaterRefFlag": "N/A",
                "deNovoAccAltGreaterClosestRefFlag": "N/A",
                "deNovoAccAltGreaterClosestAltFlag": "N/A",
                "deNovoAccFrameshiftFlag": "N/A",
                "spliceSite": 0,
                "spliceRescue": "N/A",
                "spliceFlag": deNovoDonorData["spliceFlag"],
                "frameshiftFlag": "N/A",
                "inExonicPortionFlag": "N/A",
                "CIDomainInRegionFlag": "N/A",
                "isDivisibleFlag": "N/A",
                "lowMESFlag": "N/A"}


def getPriorProbInUTRSNS(variant, boundaries, genome, transcript):
    """
    Given a variant and boundaries (either "priors" or "enigma"),
    Checks that variant is a SNS variant in a UTR
    Determines prior prob based on location (5'/3' UTR and intron/exon)
    Genome is a SequenceFileDB for genome and transcript is a pyhgvs transcript object)
       both genome and transcript are necessary to convert from genomic to transcript coordinates
    Returns a dictionary containing applicable prior and predicted qualitative enigma class
    Dictionary also contains de novo donor/acceptor ref and alt scores if applicable
    AND a spliceFlag which is equal to 1 if variant creates a better de novo splice site than ref sequence
    Rest of values in dictionary are equal to 0, "-", or N/A because they are not relevant to variant
    """
    # TO DO still need to account for creation of alternate ATG codons in 5' UTRs
    varLoc = compute.getVarLocation(variant, boundaries)
    varType = extract.getVarType(variant)
    if varLoc == "UTR_variant" and varType == "substitution":
        deNovoAccData = {"priorProb": "N/A",
                         "refMaxEntScanScore": "N/A",
                         "altMaxEntScanScore": "N/A",
                         "refZScore": "N/A",
                         "altZScore": "N/A",
                         "refSeq": "N/A",
                         "altSeq": "N/A",
                         "varStart": "N/A",
                         "varLength": "N/A",
                         "exonStart": "N/A",
                         "intronStart": "N/A",
                         "genomicSplicePos": "N/A",
                         "transcriptSplicePos": "N/A",
                         "closestGenomicSplicePos": "N/A",
                         "closestTranscriptSplicePos": "N/A",
                         "closestRefMaxEntScanScore": "N/A",
                         "closestRefZScore": "N/A",
                         "closestRefSeq": "N/A",
                         "closestAltMaxEntScanScore": "N/A",
                         "closestAltZScore": "N/A",
                         "closestAltSeq": "N/A",
                         "closestExonStart": "N/A",
                         "closestIntronStart": "N/A",
                         "altGreaterRefFlag": "N/A",
                         "altGreaterClosestRefFlag": "N/A",
                         "altGreaterClosestAltFlag": "N/A",
                         "frameshiftFlag": "N/A"}
        deNovoDonorData = {"priorProb": "N/A",
                           "refMaxEntScanScore": "N/A",
                           "altMaxEntScanScore": "N/A",
                           "refZScore": "N/A",
                           "altZScore": "N/A",
                           "refSeq": "N/A",
                           "altSeq": "N/A",
                           "varStart": "N/A",
                           "varLength": "N/A",
                           "exonStart": "N/A",
                           "intronStart": "N/A",
                           "genomicSplicePos": "N/A",
                           "transcriptSplicePos": "N/A",
                           "closestGenomicSplicePos": "N/A",
                           "closestTranscriptSplicePos": "N/A",
                           "closestRefMaxEntScanScore": "N/A",
                           "closestRefZScore": "N/A",
                           "closestRefSeq": "N/A",
                           "closestAltMaxEntScanScore": "N/A",
                           "closestAltZScore": "N/A",
                           "closestAltSeq": "N/A",
                           "closestExonStart": "N/A",
                           "closestIntronStart": "N/A",
                           "altGreaterRefFlag": "N/A",
                           "altGreaterClosestRefFlag": "N/A",
                           "altGreaterClosestAltFlag": "N/A",
                           "frameshiftFlag": "N/A"}
        varCons = extract.getVarConsequences(variant)
        if "3_prime_UTR_variant" in varCons:
            applicablePrior = LOW_PROBABILITY
            applicableClass = extract.getEnigmaClass(applicablePrior)
            spliceFlag = 0
        elif "5_prime_UTR_variant" in varCons:
            if compute.varInExon(variant):
                deNovoDonorData = getPriorProbDeNovoDonorSNS(variant, boundaries, STD_EXONIC_PORTION, genome,
                                                             transcript,
                                                             deNovoDonorInRefAcc=False)
                applicablePrior = deNovoDonorData["priorProb"]
                applicableClass = deNovoDonorData["enigmaClass"]
                spliceFlag = 0
                if compute.varInSpliceRegion(variant, donor=False, deNovo=True):
                    deNovoAccData = getPriorProbDeNovoAcceptorSNS(variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH,
                                                                  genome, transcript)
        else:
            # to account for variants in 5' UTR that are classified as other variant types
            deNovoDonorData = getPriorProbIntronicDeNovoDonorSNS(variant, genome, transcript)
            spliceFlag = deNovoDonorData["spliceFlag"]
            if spliceFlag == 1:
                applicablePrior = deNovoDonorData["priorProb"]
                applicableClass = deNovoDonorData["enigmaClass"]
            else:
                applicablePrior = LOW_PROBABILITY
                applicableClass = extract.getEnigmaClass(applicablePrior)

        return {"applicablePrior": applicablePrior,
                "applicableEnigmaClass": applicableClass,
                "proteinPrior": "N/A",
                "refDonorPrior": "N/A",
                "deNovoDonorPrior": deNovoDonorData["priorProb"],
                "refRefDonorMES": "N/A",
                "refRefDonorZ": "N/A",
                "altRefDonorMES": "N/A",
                "altRefDonorZ": "N/A",
                "refRefDonorSeq": "N/A",
                "altRefDonorSeq": "N/A",
                "refDonorVarStart": "N/A",
                "refDonorVarLength": "N/A",
                "refDonorExonStart": "N/A",
                "refDonorIntronStart": "N/A",
                "refDeNovoDonorMES": deNovoDonorData["refMaxEntScanScore"],
                "refDeNovoDonorZ": deNovoDonorData["refZScore"],
                "altDeNovoDonorMES": deNovoDonorData["altMaxEntScanScore"],
                "altDeNovoDonorZ": deNovoDonorData["altZScore"],
                "refDeNovoDonorSeq": deNovoDonorData["refSeq"],
                "altDeNovoDonorSeq": deNovoDonorData["altSeq"],
                "deNovoDonorVarStart": deNovoDonorData["varStart"],
                "deNovoDonorVarLength": deNovoDonorData["varLength"],
                "deNovoDonorExonStart": deNovoDonorData["exonStart"],
                "deNovoDonorIntronStart": deNovoDonorData["intronStart"],
                "deNovoDonorGenomicSplicePos": deNovoDonorData["genomicSplicePos"],
                "deNovoDonorTranscriptSplicePos": deNovoDonorData["transcriptSplicePos"],
                "closestDonorGenomicSplicePos": deNovoDonorData["closestGenomicSplicePos"],
                "closestDonorTranscriptSplicePos": deNovoDonorData["closestTranscriptSplicePos"],
                "closestDonorRefMES": deNovoDonorData["closestRefMaxEntScanScore"],
                "closestDonorRefZ": deNovoDonorData["closestRefZScore"],
                "closestDonorRefSeq": deNovoDonorData["closestRefSeq"],
                "closestDonorAltMES": deNovoDonorData["closestAltMaxEntScanScore"],
                "closestDonorAltZ": deNovoDonorData["closestAltZScore"],
                "closestDonorAltSeq": deNovoDonorData["closestAltSeq"],
                "closestDonorExonStart": deNovoDonorData["closestExonStart"],
                "closestDonorIntronStart": deNovoDonorData["closestIntronStart"],
                "deNovoDonorAltGreaterRefFlag": deNovoDonorData["altGreaterRefFlag"],
                "deNovoDonorAltGreaterClosestRefFlag": deNovoDonorData["altGreaterClosestRefFlag"],
                "deNovoDonorAltGreaterClosestAltFlag": deNovoDonorData["altGreaterClosestAltFlag"],
                "deNovoDonorFrameshiftFlag": deNovoDonorData["frameshiftFlag"],
                "deNovoAccPrior": deNovoAccData["priorProb"],
                "refAccPrior": "N/A",
                "refRefAccMES": "N/A",
                "refRefAccZ": "N/A",
                "altRefAccMES": "N/A",
                "altRefAccZ": "N/A",
                "refRefAccSeq": "N/A",
                "altRefAccSeq": "N/A",
                "refAccVarStart": "N/A",
                "refAccVarLength": "N/A",
                "refAccExonStart": "N/A",
                "refAccIntronStart": "N/A",
                "refDeNovoAccMES": deNovoAccData["refMaxEntScanScore"],
                "refDeNovoAccZ": deNovoAccData["refZScore"],
                "altDeNovoAccMES": deNovoAccData["altMaxEntScanScore"],
                "altDeNovoAccZ": deNovoAccData["altZScore"],
                "refDeNovoAccSeq": deNovoAccData["refSeq"],
                "altDeNovoAccSeq": deNovoAccData["altSeq"],
                "deNovoAccVarStart": deNovoAccData["varStart"],
                "deNovoAccVarLength": deNovoAccData["varLength"],
                "deNovoAccExonStart": deNovoAccData["exonStart"],
                "deNovoAccIntronStart": deNovoAccData["exonStart"],
                "closestAccRefMES": deNovoAccData["closestRefMaxEntScanScore"],
                "deNovoAccGenomicSplicePos": deNovoAccData["genomicSplicePos"],
                "deNovoAccTranscriptSplicePos": deNovoAccData["transcriptSplicePos"],
                "closestAccGenomicSplicePos": deNovoAccData["closestGenomicSplicePos"],
                "closestAccTranscriptSplicePos": deNovoAccData["closestTranscriptSplicePos"],
                "closestAccRefZ": deNovoAccData["closestRefZScore"],
                "closestAccRefSeq": deNovoAccData["closestRefSeq"],
                "closestAccAltMES": deNovoAccData["closestAltMaxEntScanScore"],
                "closestAccAltZ": deNovoAccData["closestAltZScore"],
                "closestAccAltSeq": deNovoAccData["closestAltSeq"],
                "closestAccExonStart": deNovoAccData["closestExonStart"],
                "closestAccIntronStart": deNovoAccData["closestExonStart"],
                "deNovoAccAltGreaterRefFlag": deNovoAccData["altGreaterRefFlag"],
                "deNovoAccAltGreaterClosestRefFlag": deNovoAccData["altGreaterClosestRefFlag"],
                "deNovoAccAltGreaterClosestAltFlag": deNovoAccData["altGreaterClosestAltFlag"],
                "deNovoAccFrameshiftFlag": deNovoAccData["frameshiftFlag"],
                "spliceSite": 0,
                "spliceRescue": "N/A",
                "spliceFlag": spliceFlag,
                "frameshiftFlag": "N/A",
                "inExonicPortionFlag": "N/A",
                "CIDomainInRegionFlag": "N/A",
                "isDivisibleFlag": "N/A",
                "lowMESFlag": "N/A",
                "varConsequences": ",".join(varCons)}
