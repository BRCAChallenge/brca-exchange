#!/usr/bin/env python
"""
These are shared utilities and parameters for generating the BRCA Exchange track hub for the UCSC Genome Browser
"""
import argparse
import re

PATHOGENIC_VERY_STRONG_COLOR = "30,12,43"   # Dark purple
PATHOGENIC_STRONG_COLOR = "88,58,110"       # Deep purple
PATHOGENIC_MODERATE_COLOR = "138,111,158"  # Light purple
PATHOGENIC_SUPPORTING_COLOR = "204,186,217"  # Lavender, very light purple

BENIGN_VERY_STRONG_COLOR = "2,82,66"  # Dark cyan
BENIGN_STRONG_COLOR = "46,159,134"      # Deep cyan
BENIGN_MODERATE = "124,191,178"         # Light cyan
BENIGN_SUPPORTING_COLOR = "168,224,213" # Very light cyan, mint

NO_CODE_COLOR = "13,85,191"  # Royal blue

PATHOGENIC_COLOR = "156,6,6"  # Brick red
BENIGN_COLOR = "0,180,85"       # Bluish green
VUS_COLOR = "100,100,100"       # Slate gray       

def _add_urls(s, url=None):
    """ transform a list of URLs to hrefs """
    lines = []
    for part in s.split(","):
        part = part.strip()
        if part == "":
            continue
        if part.startswith("http"):
            label = part.split("/")[-1]
            if "=" in label:
                label = label.split("=")[-1]
            part = "<a href='%s'>%s</a>" % (part, label)
            lines.append(part)
        else:
            if url == None:
                lines.append(part)
            else:
                part = "<a href='%s%s'>%s</a>" % (url, part, part)
                lines.append(part)

    return ", ".join(lines)

def _get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", help="Path to built_with_change_types.tsv file",
                        default="output/release/built_with_change_types.tsv")
    parser.add_argument("-l", "--length_threshold", help="Length threshold for short/structural variants",
                        type=int, default=50)
    parser.add_argument("--output_hg19_var", help="Short variants: hg19 output BED file",
                        default="variants.hg19.bed")
    parser.add_argument("--output_hg38_var", help="Short variants: hg38 output BED file",
                        default="variants.hg38.bed")
    parser.add_argument("--output_hg19_sv", help="Structural variants: hg19 output BED file",
                        default="structural_variants.hg19.bed")
    parser.add_argument("--output_hg38_sv", help="Structural variants: hg38 output BED file",
                        default="structural_variants.hg38.bed")
    parser.add_argument("-a", "--auto-sql-file", help="Field definitions in AutoSQL format",
                        default="brcaExchange.as")
    args = parser.parse_args()
    return args

def displayString(contents):
    if contents is None:
        return("-")
    else:
        return(contents)


def pathogenicityToColor(pathogenicity):
    if re.search("pathogen", pathogenicity, re.IGNORECASE):
        color = PATHOGENIC_COLOR
    elif re.search("benign", pathogenicity, re.IGNORECASE):
        color = BENIGN_COLOR
    else:
        color = VUS_COLOR
    return(color)

        

def acmgCodeToColor(acmgCode):
    """
    Given an ACMG provisional evidence code, return the appropriate display color
    """
    if re.search("^B.*Supporting", acmgCode, re.IGNORECASE):
        color = BENIGN_SUPPORTING_COLOR
    elif re.search("^B.*Moderate", acmgCode, re.IGNORECASE):
        color = BENIGN_MODERATE_COLOR
    elif re.search("^B.*Strong", acmgCode, re.IGNORECASE):
        color = BENIGN_STRONG_COLOR
    elif re.search("^B.*VeryStrong", acmgCode, re.IGNORECASE):
        color = BENIGN_VERY_STRONG_COLOR
    elif re.search("$BA", acmgCode, re.IGNORECASE):
        color = BENIGN_VERY_STRONG_COLOR
    elif re.search("$BS", acmgCode, re.IGNORECASE):
        color = BENIGN_STRONG_COLOR
    elif re.search("$BM", acmgCode, re.IGNORECASE):
        color = BENIGN_MODERATE_COLOR
    elif re.search("$BP", acmgCode, re.IGNORECASE):
        color = BENIGN_SUPPORTIVE_COLOR
    elif re.search("^P.*Supporting", acmgCode, re.IGNORECASE):
        color = PATHOGENIC_SUPPORTING_COLOR
    elif re.search("^P.*Moderate", acmgCode, re.IGNORECASE):
        color = PATHOGENIC_MODERATE_COLOR
    elif re.search("^P.*Strong", acmgCode, re.IGNORECASE):
        color = PATHOGENIC_STRONG_COLOR
    elif re.search("^P.*VeryStrong", acmgCode, re.IGNORECASE):
        color = PATHOGENIC_VERY_STRONG_COLOR
    elif re.search("$PVS", acmgCode, re.IGNORECASE):
        color = PATHOGENIC_VERY_STRONG_COLOR
    elif re.search("$PS", acmgCode, re.IGNORECASE):
        color = PATHOGENIC_STRONG_COLOR
    elif re.search("$PM", acmgCode, re.IGNORECASE):
        color = PATHOGENIC_MODERATE_COLOR
    elif re.search("$PP", acmgCode, re.IGNORECASE):
        color = PATHOGENIC_SUPPORTIVE_COLOR
    else:
        color = NO_CODE_COLOR
    return(color)
        
