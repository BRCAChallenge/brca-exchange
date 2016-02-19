Source: http://evs.gs.washington.edu/EVS/

Commands for exporting EA_MAF and AA_MAF from ESP for the BRCA regions:

Uses BRCA1 GRCh37 gene bounds of chr17:41191488:41322420 and BRCA2 GRCh37 
gene bounds of chr13:32889080-32973809


espExtract.py $CGLBRCA/phase1/data/esp/ESP6500SI-V2-SSA137.GRCh38-liftover.chr17.snps_indels.vcf --start 41191488 --end 41322420 --ancestry EA > esp_ea_af.txt
espExtract.py $CGLBRCA/phase1/data/esp/ESP6500SI-V2-SSA137.GRCh38-liftover.chr13.snps_indels.vcf --start 32889080 --end 32973809 --ancestry EA >> esp_ea_af.txt

espExtract.py $CGLBRCA/phase1/data/esp/ESP6500SI-V2-SSA137.GRCh38-liftover.chr17.snps_indels.vcf --start 41191488 --end 41322420 --ancestry AA > esp_aa_af.txt
espExtract.py $CGLBRCA/phase1/data/esp/ESP6500SI-V2-SSA137.GRCh38-liftover.chr13.snps_indels.vcf --start 32889080 --end 32973809 --ancestry AA >> esp_aa_af.txt

Commands to get all the VCF fields for the BRCA1 GRCh38 gene bounds of BRCA1
of chr17:43044295-43125483 and BRCA2 of chr13:32315473-32400266

espExtract.py ESP6500SI-V2-SSA137.GRCh38-liftover.chr17.snps_indels.vcf --start 43044295 --end 43125483 --full 1 > esp.brca.vcf
espExtract.py ESP6500SI-V2-SSA137.GRCh38-liftover.chr13.snps_indels.vcf --start 32315473 --end 32400266 --full 1 >> esp.brca.vcf
