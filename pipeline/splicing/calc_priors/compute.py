import calcMaxEntScanMeanStd
from calc_priors.constants import STD_ACC_INTRONIC_LENGTH, STD_ACC_EXONIC_LENGTH, STD_DE_NOVO_LENGTH, \
    STD_DONOR_INTRONIC_LENGTH, STD_DONOR_EXONIC_LENGTH, brca1CIDomains, brca2CIDomains, STD_EXONIC_PORTION, \
    STD_DONOR_SIZE, STD_ACC_SIZE
from calc_priors import extract
from calc_priors import verify


def varInExon(variant):
    """
    Given a variant, determines if variant genomic position is inside transcript boundaries
    AND if variant is in an exon
    Returns true if variant is in an exon
    """
    varOutBounds = verify.varOutsideBoundaries(variant)
    if not varOutBounds:
        varGenPos = int(variant["Pos"])
        varExons = extract.getExonBoundaries(variant)
        varStrand = verify.getVarStrand(variant)
        for exon in varExons.keys():
            exonStart = int(varExons[exon]["exonStart"])
            exonEnd = int(varExons[exon]["exonEnd"])
            if varStrand == "+":
                if varGenPos > exonStart and varGenPos <= exonEnd:
                    return True
            else:
                if varGenPos <= exonStart and varGenPos > exonEnd:
                    return True
    return False


def varInSpliceRegion(variant, donor=False, deNovo=False):
    """
    Given a variant, determines if a variant is in reference transcript's splice donor/acceptor region
    If donor=True and deNovo=False, checks if variant is in a reference splice donor region
    If donor=True and deNovo=True, checks if variant is in a de novo splice donor region
    If donor=False and deNovo=False, checks if variant is in a reference splice acceptor region
    If donor=False and deNovo=True, checks if variant is in a de novo splice acceptor region
    Returns True if variant is in a splice region, false otherwise
    """
    if donor == False and deNovo == False:
        regionBounds = extract.getSpliceAcceptorBoundaries(variant, STD_ACC_INTRONIC_LENGTH, STD_ACC_EXONIC_LENGTH)
    elif donor == False and deNovo == True:
        regionBounds = extract.getSpliceAcceptorBoundaries(variant, STD_ACC_INTRONIC_LENGTH, STD_DE_NOVO_LENGTH)
    elif donor:
        # gets reference donor splice boundaries, if deNovo = True then entireity of exon will be included below
        regionBounds = extract.getRefSpliceDonorBoundaries(variant, STD_DONOR_INTRONIC_LENGTH, STD_DONOR_EXONIC_LENGTH)
    for exon in regionBounds.keys():
        if not donor:
            regionStart = regionBounds[exon]["acceptorStart"]
            regionEnd = regionBounds[exon]["acceptorEnd"]
        else:
            regionStart = regionBounds[exon]["donorStart"]
            regionEnd = regionBounds[exon]["donorEnd"]
        withinBoundaries = verify.checkWithinBoundaries(verify.getVarStrand(variant), int(variant["Pos"]), regionStart, regionEnd)
        if withinBoundaries == True and donor == False:
            return True
        elif donor == True and deNovo == False and withinBoundaries == True:
            return True
        # because de novo donor region includes reference splice donor region and entirity of exon
        elif donor == True and deNovo == True and (withinBoundaries == True or varInExon(variant) == True):
            return True
    return False


def varInExonicPortion(variant, exonicPortionSize, deNovoLength, donor=True, deNovoDonorInRefAcc=False):
    """
    Given a variant, determines if variant in in the exonic portion as specified
    exonicPortionLength refers to the number of bases that are considered to be in the exon
    deNovoLength refers to the number of bases in the exon that are considered part of deNovo acceptor region
    if donor=True and exonicPortionSize=3, determines if variant is in first 3 bp of highest scoring window
    if donor=False and exonicPortionSize=3, determines if variant is in last 3 bp of highest scoring window
    If deNovoDonorInRefAcc=True, function is used in context of looking for de novo donor scores in ref splice acceptor sites
    if deNovoDonorInRefAcc=False, function is not used for de novo donors in ref splice acceptor sites
    Returns true if variant is in exonic portion, False otherwise
    """
    slidingWindowInfo = getMaxMaxEntScanScoreSlidingWindowSNS(variant, exonicPortionSize=exonicPortionSize,
                                                              deNovoLength=deNovoLength, donor=donor,
                                                              deNovoDonorInRefAcc=deNovoDonorInRefAcc)
    if slidingWindowInfo["inExonicPortion"]:
        return True
    return False


def isSplicingWindowInFrame(variant, exonicPortionSize, intronicPortionSize, deNovoDonorInRefAcc=False, donor=True):
    """
    Given a variant, determines ref and alt exon length and compares them
    exonicPortionSize refers to length in bp that is considered to be in exonic portion of splice site
    intronicPortionSize referes to length in bp that is considered to be in intronic portion of splice site
    deNovoDonorInRefAcc=True if looking for deNovoDonor in ref acceptor site, False otherwise
    Donor argument determines if function is used for de novo donor (donor=True) or de novo acceptor (donor=False)
    If ref and alt exon are in the same reading frame, returns True
    """
    refLength = getRefExonLength(variant, donor=donor)
    altLength = getAltExonLength(variant, exonicPortionSize, intronicPortionSize,
                                 deNovoDonorInRefAcc=deNovoDonorInRefAcc, donor=donor)
    return verify.compareRefAltExonLengths(refLength, altLength)


def isDeNovoWildTypeSplicePosDistanceDivisibleByThree(variant, exonicPortionSize, intronicPortionSize,
                                                      deNovoDonorInRefAcc=False, donor=True):
    """
    Given a variant, compares de novo splicing position with wild-type splicing position
    exonicPortionSize refers to length in bp that is considered to be in exonic portion of splice site
    intronicPortionSize referes to length in bp that is considered to be in intronic portion of splice site
    deNovoDonorInRefAcc argument=True if looking for de novo donor in reference splice acceptor region, False otherwise
    Donor argument determines if function is used for de novo donor (donor=True) or de novo acceptor (donor=False)
    If distance between de novo and wild-type donors is divisible by 3, returns True
    returns False otherwise
    This function is another way to check if a de novo splice site would cause a frameshift mutation
       If it returns True, then de novo splicing would not cause a frameshift
       If it returns False, then de novo splicing would cause a frameshift
    """
    if varInExon(variant):
        varExonNum = getVarExonNumberSNS(variant)
    else:
        if varInSpliceRegion(variant, donor=donor, deNovo=False):
            spliceBounds = getVarSpliceRegionBounds(variant, donor=donor, deNovo=False)
            varExonNum = spliceBounds["exonName"]
        else:
            varExonNum = getClosestExonNumberIntronicSNS(variant, "enigma", donor=donor)
    varStrand = verify.getVarStrand(variant)
    refExonBounds = extract.getExonBoundaries(variant)
    slidingWindowInfo = getMaxMaxEntScanScoreSlidingWindowSNS(variant, exonicPortionSize, STD_DE_NOVO_LENGTH,
                                                              donor=donor,
                                                              deNovoDonorInRefAcc=deNovoDonorInRefAcc)
    deNovoSplicePos = extract.getNewSplicePosition(variant["Pos"], varStrand, slidingWindowInfo["varWindowPosition"],
                                           slidingWindowInfo["inExonicPortion"], exonicPortionSize, intronicPortionSize,
                                           donor=donor)
    if donor:
        wildTypeSplicePos = refExonBounds[varExonNum]["exonEnd"]
        if varStrand == "+":
            distanceBetween = wildTypeSplicePos - deNovoSplicePos
        else:
            # +1 for minus strand donor because splice donor position is to the right of splice cut position
            distanceBetween = deNovoSplicePos - (wildTypeSplicePos + 1)
    else:
        wildTypeSplicePos = refExonBounds[varExonNum]["exonStart"]
        if varStrand == "+":
            distanceBetween = abs(deNovoSplicePos - wildTypeSplicePos)
        else:
            # +1 for minus strand acceptor because splice acceptor position is to the left of splice cut position
            distanceBetween = abs((wildTypeSplicePos + 1) - deNovoSplicePos)

    if distanceBetween % 3 == 0:
        return True
    return False


def varInIneligibleDeNovoExon(variant, donor=True):
    """
    Given a variant and donor argument:
        (donor=True, if for a de novo donor, donor=False, if for a de novo acceptor)
    Determines whether that variant is in an eligible exon to be evaluated for de novo splicing
    If in ineligible exon, returns True, returns False otherwise
    """
    if varInExon(variant):
        varGene = variant["Gene_Symbol"]
        varExon = getVarExonNumberSNS(variant)
        if donor:
            # last exon not eligible for de novo splice donor
            if varGene == "BRCA1" and varExon == "exon24":
                return True
            elif varGene == "BRCA2" and varExon == "exon27":
                return True
        else:
            # first exon not eligible for de novo splice acceptor
            if varExon == "exon1":
                return True
        return False


def varInCIDomain(variant, boundaries):
    """
    Given a variant, determines if variant is in a clinically important domain
    Second argument determiens which boundaries (ENIGMA or PRIORS) are used for CI domains
    Returns True if variant in CI domain
    """
    varGenPos = int(variant["Pos"])
    varGene = variant["Gene_Symbol"]
    varStrand = verify.getVarStrand(variant)
    inExon = varInExon(variant)
    if inExon:
        if varGene == "BRCA1":
            for domain in brca1CIDomains[boundaries].keys():
                domainStart = brca1CIDomains[boundaries][domain]["domStart"]
                domainEnd = brca1CIDomains[boundaries][domain]["domEnd"]
                withinBoundaries = verify.checkWithinBoundaries(varStrand, varGenPos, domainStart, domainEnd)
                if withinBoundaries:
                    return True
        elif varGene == "BRCA2":
            for domain in brca2CIDomains[boundaries].keys():
                domainStart = brca2CIDomains[boundaries][domain]["domStart"]
                domainEnd = brca2CIDomains[boundaries][domain]["domEnd"]
                withinBoundaries = verify.checkWithinBoundaries(varStrand, varGenPos, domainStart, domainEnd)
                if withinBoundaries:
                    return True
    return False


def getVarExonNumberSNS(variant):
    """
    Given a SNS variant, checks that variant is in an exon
    If variant in an exon, returns the number of the exon variant is located within in format "exonN"
    """
    if varInExon(variant):
        varGenPos = int(variant["Pos"])
        varExons = extract.getExonBoundaries(variant)
        varStrand = verify.getVarStrand(variant)
        for exon in varExons.keys():
            exonStart = varExons[exon]["exonStart"]
            exonEnd = varExons[exon]["exonEnd"]
            if varStrand == "+":
                if varGenPos > exonStart and varGenPos <= exonEnd:
                    return exon
            else:
                if varGenPos <= exonStart and varGenPos > exonEnd:
                    return exon


def getVarLocation(variant, boundaries):
    """
    Given a variant, returns the variant location as below
    Second argument is for CI domain boundaries (PRIORS or ENIGMA)
    """
    varOutBounds = verify.varOutsideBoundaries(variant)
    if varOutBounds:
        return "outside_transcript_boundaries_variant"
    inExon = varInExon(variant)
    inSpliceDonor = varInSpliceRegion(variant, donor=True, deNovo=False)
    inSpliceAcceptor = varInSpliceRegion(variant, donor=False, deNovo=False)
    if inExon:
        inCIDomain = varInCIDomain(variant, boundaries)
        if inCIDomain == True and inSpliceDonor == True:
            return "CI_splice_donor_variant"
        if inCIDomain == True and inSpliceAcceptor == True:
            return "CI_splice_acceptor_variant"
        if inCIDomain:
            return "CI_domain_variant"
        if inSpliceDonor:
            return "splice_donor_variant"
        if inSpliceAcceptor:
            return "splice_acceptor_variant"
        inGreyZone = verify.varInGreyZone(variant)
        if inGreyZone:
            return "grey_zone_variant"
        afterGreyZone = verify.varAfterGreyZone(variant)
        if afterGreyZone:
            return "after_grey_zone_variant"
        inUTR = verify.varInUTR(variant)
        if inUTR:
            return "UTR_variant"
        return "exon_variant"
    else:
        if inSpliceDonor:
            return "splice_donor_variant"
        if inSpliceAcceptor:
            return "splice_acceptor_variant"
        inUTR = verify.varInUTR(variant)
        if inUTR:
            return "UTR_variant"
        return "intron_variant"


def getVarSpliceRegionBounds(variant, donor=False, deNovo=False):
    """
    Given a variant, checks if variant is in a splice donor/acceptor region
    If donor=True, checks if variant is in a splice donor region and returns boundaries for splice donor region
      *function CANNOT be used to return de novo donor splice region bounds*
    If donor=False and deNovo=False, checks if variant is in a ref splice acceptor region and returns boundaries for splice acceptor region
    If donor=False and deNovo=True, checks if variant is in a de novo splice acceptor region and returns boundaries for that region
    If variant is in a splice region, returns a dictionary with region boundaries where variant is located
    """
    if varInSpliceRegion(variant, donor=donor, deNovo=deNovo):
        if not donor:
            if not deNovo:
                regionBounds = extract.getSpliceAcceptorBoundaries(variant, STD_ACC_INTRONIC_LENGTH, STD_ACC_EXONIC_LENGTH)
            else:
                regionBounds = extract.getSpliceAcceptorBoundaries(variant, STD_ACC_INTRONIC_LENGTH, STD_DE_NOVO_LENGTH)
            regionStartKey = "acceptorStart"
            regionEndKey = "acceptorEnd"
        else:
            regionBounds = extract.getRefSpliceDonorBoundaries(variant, STD_DONOR_INTRONIC_LENGTH, STD_DONOR_EXONIC_LENGTH)
            regionStartKey = "donorStart"
            regionEndKey = "donorEnd"
        for exon in regionBounds.keys():
            regionStart = regionBounds[exon][regionStartKey]
            regionEnd = regionBounds[exon][regionEndKey]
            withinBoundaries = verify.checkWithinBoundaries(verify.getVarStrand(variant), int(variant["Pos"]), regionStart, regionEnd)
            if withinBoundaries:
                return {"exonName": exon,
                        regionStartKey: regionStart,
                        regionEndKey: regionEnd}


def getDeNovoFrameshiftAndCIStatus(variant, boundaries, donor=True, deNovoDonorInRefAcc=False):
    """
    Given a variant, boundaries (enigma or priors), donor argument, and deNovoDonorInRefAcc argument:
      donor argument = True for de novo donors, False for de novo acceptors
      deNovoDonorInRefAcc argument = True if lookign for de novo donor in ref acceptor site, False otherwise
    Determines if new splice position causes a frameshift and would disrupt a CI Domain
    If de novo splicing would cause a frameshift, returns False
    Else, checks to see if new splice position would splice out (skip) a CI domain
    If variant de novo splice position does not cause a frameshift and does not disrupt a CI domain, reutrns True
      Returns False otherwise
    """
    frameshiftStatus = getDeNovoSpliceFrameshiftStatus(variant, donor=donor, deNovoDonorInRefAcc=deNovoDonorInRefAcc)
    # checks to make sure that variant does not cause a frameshift
    if frameshiftStatus:
        return False
    else:
        # determine if CI domain is in region that would be skipped by new splicing
        if varInExon(variant):
            varExonNum = getVarExonNumberSNS(variant)
        else:
            if varInSpliceRegion(variant, donor=donor, deNovo=False):
                spliceBounds = getVarSpliceRegionBounds(variant, donor=donor, deNovo=False)
                varExonNum = spliceBounds["exonName"]
            else:
                if donor:
                    # if a variant is in an intron de novo donor cannot splice out any of the exon
                    # so no part of a CI domain will be spliced out
                    return True
        # varExonNum is a string in the format "exonN"
        varWindowPos = getVarWindowPosition(variant, donor=donor, deNovoDonorInRefAcc=deNovoDonorInRefAcc)
        inExonicPortion = varInExonicPortion(variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH, donor=donor,
                                             deNovoDonorInRefAcc=deNovoDonorInRefAcc)
        regionStart = extract.getNewSplicePosition(variant["Pos"], verify.getVarStrand(variant), varWindowPos, inExonicPortion,
                                           STD_EXONIC_PORTION, STD_ACC_INTRONIC_LENGTH, donor=donor)
        if donor:
            # nextExonNum parses out N from varExonNum and adds 1 to get next exon number "exonN+1"
            # uses [4:] to remove "exon" from "exonN" so can add 1 to N to get N+1
            nextExonNum = "exon" + str(int(varExonNum[4:]) + 1)
            # skips to exon 5 for any variants in BRCA1 exon 3 because exon 4 does not exist in BRCA1 RefSeq transcript
            if variant["Gene_Symbol"] == "BRCA1" and nextExonNum == "exon4":
                nextExonNum = "exon5"
            refSpliceAccBounds = extract.getSpliceAcceptorBoundaries(variant, STD_ACC_INTRONIC_LENGTH, STD_ACC_EXONIC_LENGTH)
            regionEnd = refSpliceAccBounds[nextExonNum]["acceptorStart"]
        else:
            # prevExonNum parses out N from varExonNum and adds 1 to get previous exon number "exonN-1"
            # uses [4:] to remove "exon" from "exonN" so can subtract 1 to N to get N-1
            prevExonNum = "exon" + str(int(varExonNum[4:]) - 1)
            if variant["Gene_Symbol"] == "BRCA1" and prevExonNum == "exon4":
                prevExonNum = "exon3"
            refSpliceDonorBounds = extract.getRefSpliceDonorBoundaries(variant, STD_DONOR_INTRONIC_LENGTH,
                                                               STD_DONOR_EXONIC_LENGTH)
            regionEnd = refSpliceDonorBounds[prevExonNum]["donorEnd"]
        CIDomainInRegion = verify.isCIDomainInRegion(regionStart, regionEnd, boundaries, variant["Gene_Symbol"])
        if not CIDomainInRegion:
            return True
        return False


def getMaxMaxEntScanScoreSlidingWindowSNS(variant, exonicPortionSize, deNovoLength, donor=True, deNovo=False,
                                          deNovoDonorInRefAcc=False):
    """
    Given a variant, determines the maximum alt MaxEntScan score in
       a sliding window of size STD_DONOR_SIZE with the variant in each position (1-STD_DONOR_SIZE) if donor = True
       a sliding window of size STD_ACC_SIZE with the variant in each position (1-STD_ACC_SIZE) if donor = False
    This function should be used to determine window in which de novo splicing is most likely to occur
    Function can only return highest scoring window details for either de novo donor OR de novo acceptor, not both
    If donor=True, function determines highest scoring window for potential de novo donor
    If donor=False, function determines highest scoring window for potential de novo acceptor
    Returns a dictionary containing the ref and alt MaxEntScan score and z-score and position of variant for the highest scoring window
    Ref and alt seqs for the highest scoring window are included in dictionary along with varStart (0-based index of variant for formatting)
       and varLength (equal to 1 for this function, becuase this function only works for single nucleotide substitution variants
    Dictionary also containing value "inExonicPortion" that has value either True or False
       If inExonicPortion = True, then variant is in length of bp specified by exonicPortionSize of highest scoring sliding window
       If inExonicPortion = False, then variant is NOT in length of bp specified by exonicPortionSize highest scoring sliding window
    deNovoLength refers to the length of the exonic portion of a de novo splice acceptor
    deNovoDonorInRefAcc = False if NOT checking for de novo splice donors in reference splice acceptor sites
    deNovoDonorInRefAcc = True if checking for de novo splice donors in reference splice acceptor sites
    """
    if donor:
        # uses default window size for a splice donor region
        slidingWindowInfo = extract.getMaxEntScanScoresSlidingWindowSNS(variant, STD_DONOR_SIZE, donor=donor)
    else:
        # uses default window size for a splice acceptor region
        slidingWindowInfo = extract.getMaxEntScanScoresSlidingWindowSNS(variant, STD_ACC_SIZE, donor=donor)

    windowAltMaxEntScanScores = slidingWindowInfo["windowAltMaxEntScanScores"]
    # checks to see if variatn is within reference splice donor region

    inRefSpliceDonorRegion = varInSpliceRegion(variant, donor=True, deNovo=False)
    # checks to see if variant is within reference splice acceptor region
    inRefSpliceAccRegion = varInSpliceRegion(variant, donor=False, deNovo=False)
    # if variant in ref splice donor region (for de novo donor) or in ref splice acceptor region (for de novo acceptor),
    # then need to remove native splicing window from consideration for highest scoring window
    if (inRefSpliceDonorRegion == True or inRefSpliceAccRegion == True) and deNovoDonorInRefAcc == False:
        if donor:
            refSpliceBounds = getVarSpliceRegionBounds(variant, donor=donor, deNovo=False)
            if verify.getVarStrand(variant) == "+":
                refSpliceSeq = extract.getFastaSeq(extract.getVarChrom(variant), refSpliceBounds["donorStart"],
                                           refSpliceBounds["donorEnd"], plusStrandSeq=True)
            else:
                refSpliceSeq = extract.getFastaSeq(extract.getVarChrom(variant), refSpliceBounds["donorStart"],
                                           refSpliceBounds["donorEnd"], plusStrandSeq=False)
        else:
            refSpliceBounds = getVarSpliceRegionBounds(variant, donor=donor, deNovo=True)
            deNovoOffset = deNovoLength - exonicPortionSize
            # acceptorEnd +- deNovoOffset because deNovo splice acceptor region is deNovoOffset bp longer than reference splice acceptor region
            if verify.getVarStrand(variant) == "+":
                # acceptorEnd - deNovoOffset because genomic position increases from left to right on plus strand, refSeq reduced to correct length
                refSpliceSeq = extract.getFastaSeq(extract.getVarChrom(variant), refSpliceBounds["acceptorStart"],
                                           (refSpliceBounds["acceptorEnd"] - deNovoOffset), plusStrandSeq=True)
            else:
                # acceptorEnd + deNovoOffset because genomic position decreases from left to right on minus strand, refSeq reduced to correct length
                refSpliceSeq = extract.getFastaSeq(extract.getVarChrom(variant), refSpliceBounds["acceptorStart"],
                                           (refSpliceBounds["acceptorEnd"] + deNovoOffset), plusStrandSeq=False)
        for position, seqs in slidingWindowInfo["windowSeqs"].iteritems():
            if seqs["refSeq"] == refSpliceSeq:
                refSpliceWindow = position
                # removes reference splice window so it is not considered for de novo splicing
                del windowAltMaxEntScanScores[refSpliceWindow]
    # to get tuple containing sequence with variant position with maximum alt MaxEntScan score
    maxAltWindowScore = max(windowAltMaxEntScanScores.items(), key=lambda k: k[1])
    maxVarPosition = maxAltWindowScore[0]
    maxScores = slidingWindowInfo["windowScores"][maxVarPosition]
    maxSeqs = slidingWindowInfo["windowSeqs"][maxVarPosition]

    # determines if variant is in the exonic portion specified by exonicPortionLength
    inExonicPortion = False
    if donor:
        # determines if variant is in first exonicPortionSize bp of the donor region
        if maxVarPosition <= exonicPortionSize:
            inExonicPortion = True
    else:
        # determines if variant is in the last exonicPortionSize bp of the acceptor region
        if (STD_ACC_SIZE - maxVarPosition) < exonicPortionSize:
            inExonicPortion = True

    return {"refMaxEntScanScore": maxScores["refMaxEntScanScore"],
            "refZScore": maxScores["refZScore"],
            "altMaxEntScanScore": maxScores["altMaxEntScanScore"],
            "altZScore": maxScores["altZScore"],
            "refSeq": maxSeqs["refSeq"],
            "altSeq": maxSeqs["altSeq"],
            "varStart": maxVarPosition - 1,
            "varLength": 1,
            "varWindowPosition": maxVarPosition,
            "inExonicPortion": inExonicPortion}


def getClosestExonNumberIntronicSNS(variant, boundaries, donor=True):
    """
    Given a variant and boundaries (either priors or enigma),
    1. Checks that variant is in an intron or UTR and is a SNS variant
    2. Determines the exon end that is closest to that variant
    Returns the closest exon end in the format "exonN"
    If variant is not in an intron or UTR, returns "exon0"
    """
    varLoc = getVarLocation(variant, boundaries)
    if (varLoc == "intron_variant" or varLoc == "UTR_variant") and extract.getVarType(variant) == "substitution" and varInExon(
            variant) == False:
        exonBounds = extract.getExonBoundaries(variant)
        varGenPos = variant["Pos"]
        exonIntronDiffs = {}
        for exon in exonBounds.keys():
            if verify.getVarStrand(variant) == "+":
                if donor:
                    exonIntronDiff = int(varGenPos) - int(exonBounds[exon]["exonEnd"])
                else:
                    exonIntronDiff = int(exonBounds[exon]["exonStart"]) - int(varGenPos)
                if exonIntronDiff > 0:
                    exonIntronDiffs[exon] = exonIntronDiff
            else:
                if donor:
                    exonIntronDiff = int(exonBounds[exon]["exonEnd"]) - int(varGenPos)
                else:
                    exonIntronDiff = int(varGenPos) - int(exonBounds[exon]["exonStart"])
                if exonIntronDiff > 0:
                    exonIntronDiffs[exon] = exonIntronDiff
        closestExonInfo = min(exonIntronDiffs.items(), key=lambda k: k[1])
        return closestExonInfo[0]
    return "exon0"


def getClosestSpliceSiteScores(variant, deNovoOffset, donor=True, deNovo=False, deNovoDonorInRefAcc=False,
                               testMode=False):
    """
    Given a variant, determines scores for closest reference splice sequence
    Also returns sequence of closest reference splice site and genomic position of splice site
    deNovoOffset refers to difference between de novo acceptor length and exonic portion size
       If donor = True, looks for closest splice donor sequence
       If donor = False, looks for closest splice acceptor sequence
       If deNovo = True, accomodates for de novo splicing
         *Note only use argument deNovo=True in this function if donor=False
         *Function will not return correct sequence if donor=True and deNovo=True
    If exonic variant, returns a dictionary containing:
       MaxEntScan score, z-score, and splice site sequence for reference closest splice sequence
    If variant located in referene splice site, returns a dictionary containing:
       MaxEntScan score, z-score, and splice site sequence for that reference splice site sequence
    If intronic variant or variant in UTR, returns a dictionary containg:
       MaxEntScan score, z-score, and splice site sequence for reference closest splice site
       *Note if looking for closest ref acceptor for a variant in an intron, use deNovoOffset=0
    Return dictionary also contains necessary formatting variables for splice site sequence (exonStart, intronStart)
    deNovoDonorInRefAcc = False if NOT checking for de novo splice donor sites in reference splice acceptor sites
    deNovoDonorInRefAcc = True if checking for de novo splice donor sites in reference splice acceptor sites
    """
    varGenPos = int(variant["Pos"])
    varChrom = extract.getVarChrom(variant)
    varLoc = getVarLocation(variant, "enigma")

    if (varInExon(variant) == True and deNovo == False) or (varLoc == "intron_variant" or varLoc == "UTR_variant"):
        if varInExon(variant):
            exonNumber = getVarExonNumberSNS(variant)
            exonName = exonNumber
        if (varLoc == "intron_variant" or varLoc == "UTR_variant") and varInExon(variant) == False:
            exonNumber = getClosestExonNumberIntronicSNS(variant, "enigma", donor=donor)
            exonName = exonNumber
        if donor:
            refSpliceDonorBounds = extract.getRefSpliceDonorBoundaries(variant, STD_DONOR_INTRONIC_LENGTH,
                                                               STD_DONOR_EXONIC_LENGTH)
            closestSpliceBounds = refSpliceDonorBounds[exonNumber]
        else:
            refSpliceAccBounds = extract.getSpliceAcceptorBoundaries(variant, STD_ACC_INTRONIC_LENGTH, STD_ACC_EXONIC_LENGTH)
            closestSpliceBounds = refSpliceAccBounds[exonNumber]
    if varInSpliceRegion(variant, donor=donor, deNovo=deNovo) == True and deNovoDonorInRefAcc == False:
        closestSpliceBounds = getVarSpliceRegionBounds(variant, donor=donor, deNovo=deNovo)
        exonName = closestSpliceBounds["exonName"]
    if donor:
        if verify.getVarStrand(variant) == "+":
            refSeq = extract.getFastaSeq(varChrom, closestSpliceBounds["donorStart"], closestSpliceBounds["donorEnd"],
                                 plusStrandSeq=True)
            # splice site is 3 bp to the right of donor Start (+3 because plus strand numbering increases from left to right)
            # splice site is 3 bp to the right because exon end is 3 bp to the right of donor start
            genomicSplicePos = closestSpliceBounds["donorStart"] + 3
        else:
            refSeq = extract.getFastaSeq(varChrom, closestSpliceBounds["donorStart"], closestSpliceBounds["donorEnd"],
                                 plusStrandSeq=False)
            # splice site is 3 bp to the right of donor Start (-3 because minus strand numbering decreases from left to right)
            # splice site is 3 bp to the right because exon end is 3 bp to the right of donor start
            genomicSplicePos = closestSpliceBounds["donorStart"] - 3
        exonStart = 0
        intronStart = STD_EXONIC_PORTION
    else:
        # acceptorEnd +- deNovoOffset because deNovo splice acceptor region is deNovoOffset bp longer than reference splice acceptor region
        # for plus strand it is acceptorEnd - deNovoOffset because
        # the genomic position increases from left to right on the plus strand and subtraction reduces the refSeq to correct length
        # for minus strand it is acceptorEnd + deNovoOffset because
        # the genomic position decreases from left to right on the minus strand and addition reduces the refSeq to correct length
        if verify.getVarStrand(variant) == "+":
            refSeq = extract.getFastaSeq(varChrom, closestSpliceBounds["acceptorStart"],
                                 (closestSpliceBounds["acceptorEnd"] - deNovoOffset), plusStrandSeq=True)
            # splice site is 3 bp to the left of reference acceptor End (-3 because plus strand numbering increases from left to right)
            # minus deNovoOffset because deNovo splice acceptor region is deNovoOffset bp longer than reference splice acceptor region
            genomicSplicePos = closestSpliceBounds["acceptorEnd"] - 3 - deNovoOffset
        else:
            refSeq = extract.getFastaSeq(varChrom, closestSpliceBounds["acceptorStart"],
                                 (closestSpliceBounds["acceptorEnd"] + deNovoOffset), plusStrandSeq=False)
            # splice site is 3 bp to the left of reference acceptor End (+3 because minus strand numbering decreases from left to right)
            # plus deNovoOffset because deNovo splice acceptor region is deNovoOffset bp longer than reference splice acceptor region
            genomicSplicePos = closestSpliceBounds["acceptorEnd"] + 3 + deNovoOffset
        exonStart = len(refSeq) - STD_EXONIC_PORTION
        intronStart = 0
    if not testMode:
        # to prevent issue with running max ent scan score on unittests
        if exonName == "exon0":
            return {"exonName": "N/A",
                    "sequence": "N/A",
                    "exonStart": "N/A",
                    "intronStart": "N/A",
                    "maxEntScanScore": "N/A",
                    "zScore": "N/A",
                    "genomicSplicePos": "N/A"}
        closestMaxEntScanScore = calcMaxEntScanMeanStd.runMaxEntScan(refSeq, donor=donor)
        closestZScore = extract.getZScore(closestMaxEntScanScore, donor=donor)
        return {"exonName": exonName,
                "sequence": refSeq.upper(),
                "exonStart": exonStart,
                "intronStart": intronStart,
                "maxEntScanScore": closestMaxEntScanScore,
                "zScore": closestZScore,
                "genomicSplicePos": genomicSplicePos}
    else:
        return {"exonName": exonName,
                "sequence": refSeq.upper(),
                "genomicSplicePos": genomicSplicePos}


def getRefExonLength(variant, donor=True):
    """
    Given a variant, returns the length of the reference exon
    If variant is in an exon, returns length of that exon
    If variant is in a reference splice region, returns length of exon in which exonic portion is included
    If variant is in intron, returns exon in which either closest splice donor or acceptor is included depending on donor argument
      If donor=True, returns exon length for previous exon
      If donor=False, returns exon length for subsequent exon
    """
    if varInExon(variant):
        varExonNum = getVarExonNumberSNS(variant)
    else:
        if varInSpliceRegion(variant, donor=donor, deNovo=False):
            spliceBounds = getVarSpliceRegionBounds(variant, donor=donor, deNovo=False)
            varExonNum = spliceBounds["exonName"]
        else:
            varExonNum = getClosestExonNumberIntronicSNS(variant, "enigma", donor=donor)
    exonBounds = extract.getExonBoundaries(variant)
    if verify.getVarStrand(variant) == "-":
        varExonStart = int(exonBounds[varExonNum]["exonStart"])
        varExonEnd = int(exonBounds[varExonNum]["exonEnd"])
        # +1 is not included in the below equation for exonLength
        # because due to RefSeq numbering varExonEnd is 1 bp too long
        # varExonEnd is first intronic base (+1 position)
        # for minus strand genes
        exonLength = varExonStart - varExonEnd
    else:
        varExonStart = int(exonBounds[varExonNum]["exonStart"])
        varExonEnd = int(exonBounds[varExonNum]["exonEnd"])
        # +1 is not included in the below equatio for exonLength
        # because due to RefSeq numbering varExonStart is 1 bp too long
        # varExonStart is last intronic base (-1 position)
        # for plus strand genes
        exonLength = varExonEnd - varExonStart
    return exonLength


def getAltExonLength(variant, exonicPortionSize, intronicPortionSize, deNovoDonorInRefAcc=False, donor=True):
    """
    Given a variant and the exonic portion size,
    returns the length of the alternate exon after splicing occurs in max MES window
    Donor argument determines if function is used for de novo donor (donor=True) or de novo acceptor (donor=False)
    deNovoDonorInRefAcc=True if looking for deNovoDonor in ref acceptor site, False otherwise
    """
    if varInExon(variant):
        varExonNum = getVarExonNumberSNS(variant)
    else:
        if varInSpliceRegion(variant, donor=donor, deNovo=False):
            spliceBounds = getVarSpliceRegionBounds(variant, donor=donor, deNovo=False)
            varExonNum = spliceBounds["exonName"]
        else:
            varExonNum = getClosestExonNumberIntronicSNS(variant, "enigma", donor=donor)
    exonBounds = extract.getExonBoundaries(variant)
    slidingWindowInfo = getMaxMaxEntScanScoreSlidingWindowSNS(variant, exonicPortionSize, STD_DE_NOVO_LENGTH,
                                                              donor=donor,
                                                              deNovoDonorInRefAcc=deNovoDonorInRefAcc)
    newSplicePos = extract.getNewSplicePosition(variant["Pos"], verify.getVarStrand(variant), slidingWindowInfo["varWindowPosition"],
                                        slidingWindowInfo["inExonicPortion"], exonicPortionSize, intronicPortionSize,
                                        donor=donor)
    if verify.getVarStrand(variant) == "-":
        if donor:
            varExonStart = int(exonBounds[varExonNum]["exonStart"])
            # +1 to account for all positions including newSplicePos
            # adding one to exonStart increases length by 1 bp because numbering decreases from left to right
            exonLength = varExonStart - newSplicePos + 1
        else:
            varExonEnd = int(exonBounds[varExonNum]["exonEnd"])
            # -1 to account for all position including newSplicePos
            # subtracting one increases length by 1 bp because numbering decreases from left to right
            exonLength = newSplicePos - varExonEnd - 1
    else:
        varExonEnd = int(exonBounds[varExonNum]["exonEnd"])
        if donor:
            refExonLength = getRefExonLength(variant, donor=donor)
            # need to compare to refExonLength because of + strand gene
            exonLength = refExonLength - (varExonEnd - newSplicePos)
        else:
            exonLength = varExonEnd - newSplicePos
    return exonLength


def getDeNovoSpliceFrameshiftStatus(variant, donor=True, deNovoDonorInRefAcc=False):
    """
    Given a variant, determiens if de novo splice site (either donor or acceptor based on donor argument)
      causes a frameshift
    Returns true if variant's de novo splice site causes a frameshift, false otherwise

    deNovoDonorInRefAcc is True if looking for de novo donor in reference splice acceptor site, False otherwise
    Frameshift is determined in 2 ways: inFrame and isDivisible
    """
    # if inFrame == False then alt and ref exons are not in the same reading frame
    inFrame = isSplicingWindowInFrame(variant, STD_EXONIC_PORTION, STD_ACC_INTRONIC_LENGTH,
                                      deNovoDonorInRefAcc=deNovoDonorInRefAcc, donor=donor)
    # if isDivisble == Flase then distance between old and new splice position is not divislbe by 3
    isDivisible = isDeNovoWildTypeSplicePosDistanceDivisibleByThree(variant, STD_EXONIC_PORTION,
                                                                    STD_ACC_INTRONIC_LENGTH,
                                                                    deNovoDonorInRefAcc=deNovoDonorInRefAcc,
                                                                    donor=donor)
    if inFrame == False or isDivisible == False:
        return True
    return False


def getVarWindowPosition(variant, donor=True, deNovoDonorInRefAcc=False):
    """
    Given a variant, determines window position for highest scoring sliding window
    donor=True if function being used for splice donor, donor=False if function being used for splice acceptor
    Returns integer 1-STD_DONOR_SIZE based on variant position in highest scoring window if donor=True
    Returns integer 1-STD_ACC_SIZE based on variant position in highest scoring window if donor=False
    deNovoDonorInRefAcc=True if looking for deNovoDonor in ref acceptor site, False otherwise
    """
    slidingWindowInfo = getMaxMaxEntScanScoreSlidingWindowSNS(variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH,
                                                              donor=donor, deNovoDonorInRefAcc=deNovoDonorInRefAcc)
    varWindowPos = slidingWindowInfo["varWindowPosition"]
    return varWindowPos