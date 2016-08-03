
XMLPATH = /hive/groups/cgl/brca
DATA = /hive/groups/cgl/brca /hive/groups/cgl/brca/InvitaeStudy/data
CLINVAR_XML = $(XMLPATH)/ClinVarFullRelease_2016-01.xml
INVITAE_VCF = data/Invitae_Nov_submission_BRCA_vcf_default

all: data/ClinVarBRCA.selectedLabsOnly.txt

data/ClinVarBRCA.selectedLabsOnly.txt: data/ClinVarBRCA.plusInvitae.forSpreadsheet.txt
	head -1 $? > $@
	grep  "Ambry Genetics" $? >> $@
	grep "Sharing Clinical Reports Project (SCRP)" $? >> $@
	grep Invitae $? |sed 's/Invitae /Invitae/'>> $@
	grep GeneDx $? >> $@
	grep "Molecular Genetics Diagnostic Laboratory,Children's Hospital of Eastern Ontario" $? >> $@
	grep "Emory Genetics Laboratory" $? >> $@
	grep Counsyl $? >> $@

data/ClinVarBRCA.plusInvitae.forSpreadsheet.txt: data/ClinVarBrca.forSpreadsheet.txt $(INVITAE_VCF)
	src/mergeInvitaeVariants.py data/ClinVarBrca.forSpreadsheet.txt $(INVITAE_VCF) > $@

data/ClinVarBrca.forSpreadsheet.txt: $(XMLPATH)/ClinVarBrca.xml
	clinVarParse.py  $? > $@



$(XMLPATH)/ClinVarBrca.xml: $(CLINVAR_XML)
	python3 clinVarBrca.py > $@
