# Luigi (WIP)

Contains code used to stitch together processes related to the pipeline.

## Instructions

To run, you must first install Luigi: `pip install luigi`.

Also, make sure that all dependencies of the pipeline directory are properly installed. Currently, the Luigi scripts herein require the `brca-resources` directory described in the `brca-exchange/pipeline` dependency instructions to be located at `brca-exchange/pipeline/brca/brca-resources`.

Username and password are required to download files from BIC. They can be found in `/hive/groups/cgl/brca/phase1/data/bic/account.txt` at UCSC.

To execute, navigate to `brca-exchange/pipeline/luigi` and run `python -m luigi --module CompileVCFFiles RunAll --u {username} --p {password} --local-scheduler`.