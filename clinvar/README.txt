#
# Parse data out of the ClinVar XML file into VCF format.
# 
# 1. Download the latest ClinVar XML file from ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/
# 
# 2. Extract the BRCA variants from the larger ClinVar XML file.  LONG RUN TIME
# While doing so, filter out the blank lines that are probably introduced by
# a magic formula that keeps the code from choking on extended ascii 
# characters.  Produces a smaller XML.
#
# clinVarBrca.py $CLINVARXML | awk '{ if (NF > 0) { print }}' > $CLINVAR_BRCA

#
# 3. Extract data from the XML file into tab-delimited format
# clinVarParse.py $CLINVAR_BRCA --assembly GRCh38 > $CLINVAR_BRCA_GRCH38_TAB

#
# 4. Convert the tab data into VCF format, using convert_tsv_to_vcf from utils
# ../utils/convert_tsv_to_vcf.py -i $CLINVAR_BRCA_GRCH38_TAB -s "#GRCH38" -o $CLINVAR_BRCA_VCF
 
