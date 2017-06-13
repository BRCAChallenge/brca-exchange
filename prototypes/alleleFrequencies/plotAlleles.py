import vcfAlleleLibrary
import alleleRenderLib
import gzip
import shutil

# with gzip.open('1kGSample.vcf', 'rb') as f_in, gzip.open('1kGSample.vcf.gz', 'wb') as f_out:
#     shutil.copyfileobj(f_in, f_out)

# with open('1kGSample.vcf', 'rb') as f_in, gzip.open('1kGSample.vcf.gz', 'wb') as f_out:
#     shutil.copyfileobj(f_in, f_out)

filename1 = '1kGSample.vcf.gz'
alleleList = vcfAlleleLibrary.vcfFetchAlleles(filename1, 13, 32316852, 32316853)
for allele in alleleList:
    print vars(allele)
    for subpop in allele.pops:
        print vars(subpop)
alleleRenderLib.plotAlleles(alleleList, filename1)

filename2 = 'testExac.vcf.gz'
alleleList = vcfAlleleLibrary.vcfFetchAlleles(filename2, 13, 32316434, 32316435)
for allele in alleleList:
    print vars(allele)
    for subpop in allele.pops:
        print vars(subpop)
alleleRenderLib.plotAlleles(alleleList, filename2)