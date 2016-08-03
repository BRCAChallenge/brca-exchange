BRCA1 = {"hg38": {"start": 43000000,
                  "sequence": open("../resources/brca1_hg38.txt", "r").read()},
         "hg19": {"start": 41100000,
                  "sequence": open("../resources/brca1_hg19.txt", "r").read()}}
BRCA2 = {"hg38": {"start": 32300000,
                  "sequence": open("../resources/brca2_hg38.txt", "r").read()},
         "hg19": {"start": 32800000,
                  "sequence": open("../resources/brca2_hg19.txt", "r").read()}}

def ref_correct(v, version="hg38"):
    chr, pos, ref, alt = v
    pos = int(pos)
    if chr == "13":
        seq = BRCA2[version]["sequence"]
        pos = pos - 1 - BRCA2[version]["start"]
    elif chr == "17":
        seq = BRCA1[version]["sequence"]
        pos = pos - 1 - BRCA1[version]["start"]
    else:
        assert(False)

    genomeRef = seq[pos:pos+len(ref)].upper()
    if len(ref) != 0 and len(genomeRef)==0:
        print v
        raise Exception("ref not inside BRCA1 or BRCA2")
    if (genomeRef != ref):
        return False
    else:
        return True

    
def variant_equal(v1, v2, version="hg38"):
    " return (edited1, edited2) "
    if v1 == v2:
        return True

    chr1, pos1, ref1, alt1 = v1
    chr2, pos2, ref2, alt2 = v2
    pos1 = int(pos1)
    pos2 = int(pos2)

    if chr1 != chr2:
        return False
    if (len(ref1) - len(alt1)) != (len(ref2) - len(alt2)):
        return False

    # if len(ref2)>100 or len(ref1)>100:
    #     return False

    # make sure that v1 is upstream of v2
    if pos1 > pos2:
        return variant_equal(v2, v1)

    # lift coordinates and make everything 0-based
    if chr1 == "13":
        seq = BRCA2[version]["sequence"]
        pos1 = pos1 - 1 - BRCA2[version]["start"]
        pos2 = pos2 - 1 - BRCA2[version]["start"]
    elif chr1 == "17":
        seq = BRCA1[version]["sequence"]
        pos1 = pos1 - 1 - BRCA1[version]["start"]
        pos2 = pos2 - 1 - BRCA1[version]["start"]
    else:
        assert(False)

    # replace vcf ref string with alt string
    edited_v1 = seq[0:pos1]+alt1+seq[pos1+len(ref1):]
    edited_v2 = seq[0:pos2]+alt2+seq[pos2+len(ref2):]

    return edited_v1 == edited_v2


