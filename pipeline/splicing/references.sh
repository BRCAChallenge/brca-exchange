#!/usr/bin/env bash
# Download and verify references - to be called from within the docker
echo "Downloading references. This may take 1-2 hours..."
wget -N -P /references http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
wget -N -P /references ftp://ftp.ensembl.org/pub/release-92/variation/VEP/homo_sapiens_vep_92_GRCh38.tar.gz
echo "Validating MD5 of references..."
md5sum -c tests/md5/references.md5
gunzip -k -f /references/hg38.fa.gz
mkdir -p /references/vep
echo "Unpacking VEP references..."
tar -xzf /references/homo_sapiens_vep_92_GRCh38.tar.gz -C /references/vep
echo "Calculating one variant to force installation of references..."
python calcVarPriors.py calc ./tests/variants_one.tsv /tmp/priors_one.tsv
echo "Reference download and installation complete."
