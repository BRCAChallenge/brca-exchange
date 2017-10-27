import os
import utils

download_data_cmd =
"""
G1K="../data"
BRCA_RESOURCES=$G1K
mkdir $G1K

echo downloading 1000 genome variant data from ftp....takes about 5 minutes....
echo -e "\ndata location: $G1K"
pushd $G1K
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
popd

echo extract BRCA gene region from chr13 and chr17
tabix -h $G1K/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 13:32889080-32973809 > $G1K/chr13_brca2_1000g_GRCh37.vcf
tabix -h $G1K/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 17:41191488-41322420 > $G1K/chr17_brca1_1000g_GRCh37.vcf
vcf-concat $G1K/chr13_brca2_1000g_GRCh37.vcf $G1K/chr17_brca1_1000g_GRCh37.vcf  > $G1K/brca12_1000g_GRCh37.vcf
CrossMap.py vcf $BRCA_RESOURCES/hg19ToHg38.over.chain.gz  $G1K/brca12_1000g_GRCh37.vcf $BRCA_RESOURCES/hg38.fa  $G1K/1000G_brca.hg38.vcf
vcf-sort $G1K/1000G_brca.hg38.vcf > $G1K/1000G_brca.sorted.hg38.vcf

echo Done!
"""

extract_high_af_cmd = 
"""
vcftools --vcf """ +
exon_fn + " --maf " + str(0.01) + """
"""

# File names
start_fn = """1000G_brca.hg38.vcf"""
exon_fn = """1000G_brca.hg38_exon.vcf"""
af_fn = """1000G_brca.hg38_exon_high_af.vcf"""

# Download data files
os.system(download_data_cmd)
# TODO RefSeq files?

# Extract exons
starts, ends = utils.exon_starts_ends('RefSeq.genepred', ['NM_000059.3', 'NM_007294.3'])
with open('exon_starts_ends.txt', 'wb') as exon_ind_file:
  for start in starts:
utils.extract_exon(start_fn, exon_fn)

# Extract high allele frequencies
os.system(extract_high_af_cmd)
