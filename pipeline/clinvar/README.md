# Parse data out of the ClinVar XML file into VCF format.

In order for the makefile to work, you will need a directory with the following structure:

- brca
  - pipeline-data
    - data
      - ClinVar
      - pipeline_input

You will also need to create an empty file called `ClinVarBrca.xml` inside the ClinVar directory `brca/pipeline-data/data/ClinVar` if you don't already have this file.

1. Make sure you have a an environment variable for the `brca/pipeline-data` directory named `BRCA_PIPELINE_DATA`. It's advisable to create this variable in your `.bashrc` file to make it readily available for future use.

2. Download the latest ClinVar XML file in gzip format from ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ (it should be titled `ClinVarFullRelease_00-latest.xml.gz`) and move it into `brca/pipeline-data/data/ClinVar`.

3. cd into `brca-exchange/pipeline/clinvar` and run `make`.

4. Find something to do while the process is running. Be patient, it takes a minute!