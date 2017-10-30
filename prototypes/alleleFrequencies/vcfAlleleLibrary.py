import vcf

class Allele(object):
    """This object an allele made from a vcf record. A novel feature of this representation is 
    that it allows subpopulation data to be extracted from the record and well organized"""
    def __init__(self, chrom, pos, ref, alt): #chrom pos ref and alt might be passed here actualyl
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.pops = []

class SubPop(object):
    """This object is a subpopulation that contains its name, allele counts, allele numbers, and allele frequency"""
    def __init__(self, namepop=None, ac=None, an=None, af=None): #for 1KG, af=#.##
        self.ac = ac
        self.an = an
        self.af = af #Can be None for ExAC
        self.namepop = namepop
        if self.ac is not None:
            self.af = float(self.ac)/float(self.an)

def vcfFetchAlleles(vcffile, chrom, start, end):
    """Converts a range of vcf file records into Allele objects"""
    vcfReader = vcf.Reader(filename=vcffile)
    alleleList = []
    #print 'hello'
    for vcfrec in vcfReader.fetch(chrom, start, end):
        alleleList.extend(vcfRecToAlleles(vcfrec))
    return alleleList

def vcfRecToAlleles(vcfrec):
    """Converts record objects to a list allele objects"""
    alleleList = []
    for altindex in xrange(len(vcfrec.ALT)):
        allele = vcfAltToAllele(vcfrec, altindex)
        alleleList.append(allele)
    return alleleList

def vcfAltToAllele(vcfrec, altindex):
    "Converts a single record to a fully characterized (?) allele object"
    allele = Allele(vcfrec.CHROM, vcfrec.POS, vcfrec.REF, str(vcfrec.ALT[altindex]))
    popList = vcfRecToPopNames(vcfrec)
    if 'NFE' in popList:
        for popname in popList:
            subpop = ExacVcfToSubpop(vcfrec, altindex, popname)
            allele.pops.append(subpop)
    else:
        for popname in popList:
            subpop = TgVcfToSubpop(vcfrec, altindex, popname)
            allele.pops.append(subpop)
    return allele

def vcfRecToPopNames(vcfrec):
    """Generates a list of population names from the record object"""
    exacPopList = []
    tgPopList = []
    info = vcfrec.INFO
    for field in info.keys():
        if field[0:3] == 'AC_' and field.isupper() and len(field) < 7:
            exacPopList.append(field[3:])
        if field[3:] == '_AF':
            tgPopList.append(field[0:3])
    if len(exacPopList) > len(tgPopList):
        popList = exacPopList
    else:
        popList = tgPopList
    return popList

def ExacVcfToSubpop(vcfrec, altindex, popname):
    """Returns a subpop object with AC and AN fields from info column of vcf"""
    acs = vcfrec.INFO["AC_" + popname]
    subpop = SubPop(namepop = popname, ac=acs[altindex], an=vcfrec.INFO["AN_" + popname]) #ACs are list type, ANs are int type
    return subpop

def TgVcfToSubpop(vcfrec, altindex, popname):
    afs =vcfrec.INFO[popname + "_AF"]
    subpop = SubPop(namepop= popname, af=afs[altindex])
    return subpop





