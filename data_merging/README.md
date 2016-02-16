Pipeline:

1. one_variant_per_row.py

take all the VCF, read through them, if there is a variant line that has mulitiple ALT field looking like this:
CHROM  POS     REF   ALT
13           12345   A        T,G
It means this row is actually two variants (it happens quite frequently), in this case, all rows like this has to be expanded into two rows:
CHROM  POS     REF   ALT
13           12345   A        T
13           12345   A        G
This one_variant_per_row.py script does the expansion and generate new vcf files for each raw vcf input

2.  repeat_merging.py
take the output from step one, merge repetitive variants within each vcf file by string comparison. Now after step 1 and 2's preprocessing of the vcf files, all the vcfs are saved as following and they are used to be fed to step 3.
ENIGMA_FILE = "../data/enigma_variants_9-29-2015.tsv"
GENOME1K_FILE = "../data/allVcf/no_repeats/1000_genomes.brca.no_sample.ovpr.no_repeats.vcf"
CLINVAR_FILE = "../data/allVcf/no_repeats/clinvar.brca.ovpr.no_repeats.vcf"
LOVD_FILE = "../data/allVcf/no_repeats/lovd.brca.ovpr.no_repeats.vcf"
EX_LOVD_FILE = "../data/allVcf/no_repeats/ex_lovd.brca.ovpr.no_repeats.vcf"
BIC_FILE = "../data/allVcf/no_repeats/bic.brca.no_repeats.vcf"
EXAC_FILE = "../data/allVcf/no_repeats/exac.brca.ovpr.no_repeats.vcf"

3 . variant_merging.py
merge the enigma file with all the vcf files 



