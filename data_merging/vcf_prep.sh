cut -f1-9 ../../Data/allVcfs/1000_genomes.brca.vcf > ../../Data/allVcfs/1000_genomes.brca.no_sample.vcf
python exac_VEP_expansion.py -i ../../Data/allVcfs/exac_BRCA12.sorted.vcf -o ../../Data/allVcfs/exac.VEP_expanded.vcf
