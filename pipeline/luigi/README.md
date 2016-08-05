# Luigi (WIP)

Contains code used to stitch together processes related to the pipeline.

## Instructions

To run, you must first install Luigi `pip install luigi`.

### Clinvar

This process downloads the latest ClinVar data and converts it to VCF format. Eventually, the script will be much more robust and require less setup.

Required directory structure is as follows:

- root
  - brca
    - pipeline-data
      - data
        - ClinVar
          - ClinVarBrca.xml
        - pipeline_input
  - brca-exchange
    - website
    - pipeline
    - deployment
    - README.md
    
Note that there must be an empty `ClinVarBrca.xml` file in the ClinVar directory to run successfully. `brca-exchange` is the repo that contains this file.

You will also need an environment variable called `BRCA_PIPELINE_DATA` that points to the `brca/pipeline-data` directory (e.g. add `export BRCA_PIPELINE_DATA="/path/to/brca/pipeline-data"` to your `.bashrc` file).

To execute, navigate to brca-exchange/pipeline/luigi and run `python -m luigi --module ConvertLatestClinvarToVCF ConvertLatestClinvarToVCF --local-scheduler`.