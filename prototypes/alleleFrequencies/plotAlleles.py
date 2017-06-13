import vcfAlleleLibrary
import alleleRenderLib


def fenderAllele(vcfName, chrom, start, end, color, title, pdfName):
    """fetches and renders alleleList"""
    alleleList = vcfAlleleLibrary.vcfFetchAlleles(vcfName, chrom, start, end)
    alleleRenderLib.plotAlleles(alleleList, color, title, pdfName)



fenderAllele('1kGSample.vcf.gz', 13, 32316852, 32316853, '1000 Genomes', 'cyan', '1000GenomesOutput.pdf')
fenderAllele('testExac.vcf.gz', 13, 32319058, 32319059, 'ExAC Alleles', 'salmon', 'ExACOutput.pdf')


