Commands for exporting EA_MAF and AA_MAF from ESP for the BRCA regions:

Uses BRCA1 gene bounds of chr17:41191488:41322420 and BRCA2 gene bounds of chr13:32889080-32973809


espExtract.py $CGLBRCA/phase1/data/esp/ESP6500SI-V2-SSA137.GRCh38-liftover.chr17.snps_indels.vcf --start 41191488 --end 41322420 --ancestry EA > esp_ea_af.txt
espExtract.py $CGLBRCA/phase1/data/esp/ESP6500SI-V2-SSA137.GRCh38-liftover.chr13.snps_indels.vcf --start 32889080 --end 32973809 --ancestry EA >> esp_ea_af.txt

espExtract.py $CGLBRCA/phase1/data/esp/ESP6500SI-V2-SSA137.GRCh38-liftover.chr17.snps_indels.vcf --start 41191488 --end 41322420 --ancestry AA > esp_aa_af.txt
espExtract.py $CGLBRCA/phase1/data/esp/ESP6500SI-V2-SSA137.GRCh38-liftover.chr13.snps_indels.vcf --start 32889080 --end 32973809 --ancestry AA >> esp_aa_af.txt

