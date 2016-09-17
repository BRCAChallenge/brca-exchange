### wrangling BRCA variant data from openSNP https://opensnp.org

#### 1. Download and unzip data
    cd data  
    wget https://opensnp.org/data/zip/opensnp_datadump.current.zip
    mkdir datadump_opensnp_9-8-2016
    unzip opensnp_datadump.current.zip -d datadump_opensnp_9_8_2016/
    
#### 2. re-organize data based on the content, (23andme, ancestry, ftdna, fitbit, picture, phenotype, readme)
    python preprocess_files.py
    
#### 3. extract BRCA region from 23andme and ancestry data while converting to build GRCh38
    python extract_brca.py

