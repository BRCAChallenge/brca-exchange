## Generate data for ESP

1. Generate a working directory and point to it with the environment variable $ESP
2. pushd $ESP; wget http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.tar.gz; popd
3. pushd $ESP; tar xzvf ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.tar.gz; popd
4. espExtract.py $ESP/ESP6500SI-V2-SSA137.GRCh38-liftover.chr17.snps_indels.vcf --start 43044295 --end 43125483 --full 1 -o $ESP/esp.brca1.vcf
5. espExtract.py $ESP/ESP6500SI-V2-SSA137.GRCh38-liftover.chr13.snps_indels.vcf --start 43044295 --end 43125483 --full 1 -o $ESP/esp.brca2.vcf
6. vcf-concat $ESP/esp.brca1.vcf $ESP/esp.brca2.vcf > $ESP/esp.brca12.hg38.vcf
7. vcf-sort $ESP/esp.brca12.hg38.vcf > $ESP/esp.brca12.sorted.hg38.vcf
8. cp $ESP/esp.brca12.sorted.hg38.vcf $PIPELINE_INPUT
