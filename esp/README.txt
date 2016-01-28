Commands for exporting EA_MAF and AA_MAF from ESP for the BRCA regions:

Uses BRCA1 gene bounds of chr17:43044295-43125483 and BRCA2 gene bounds of
chr13:32315474-32400266

espExtract.py $CGLBRCA/phase1/data/esp/ESP6500SI-V2-SSA137.GRCh38-liftover.chr17.snps_indels.vcf --start 43044295 --end 43125483 --ancestry EA > esp_ea_af.txt
espExtract.py $CGLBRCA/phase1/data/esp/ESP6500SI-V2-SSA137.GRCh38-liftover.chr13.snps_indels.vcf --start 32315474 --end 32400266 --ancestry EA >> esp_ea_af.txt

espExtract.py $CGLBRCA/phase1/data/esp/ESP6500SI-V2-SSA137.GRCh38-liftover.chr17.snps_indels.vcf --start 43044295 --end 43125483 --ancestry AA > esp_aa_af.txt
espExtract.py $CGLBRCA/phase1/data/esp/ESP6500SI-V2-SSA137.GRCh38-liftover.chr13.snps_indels.vcf --start 32315474 --end 32400266 --ancestry AA >> esp_aa_af.txt
