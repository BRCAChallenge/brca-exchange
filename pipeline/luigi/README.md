# Luigi (WIP)

Contains code used to stitch together processes related to the pipeline.

## Instructions

To run, you must first install Luigi: `pip install luigi`.

Also, make sure that all dependencies of the pipeline directory are properly installed.

Create environment variables for paths to an output directory, BRCA resources directory, and a parent directory which contains the following directories for each of the separate data sources: `BIC, ClinVar, enigma, ESP, exac, exLOVD, G1K, and LOVD.`

Username and password are required to download files from BIC. They can be found in `/hive/groups/cgl/brca/phase1/data/bic/account.txt` at UCSC.

To run: `python -m luigi --module CompileVCFFiles RunAll --u {username} --p {password} --local-scheduler --output-dir $OUTPUT_DIR --resources-dir $BRCA_RESOURCES --file-parent-dir $PARENT_DIR`

You can replace `RunAll` with individual tasks, or comment out required tasks in the `RunAll` task to control which tasks are run. A task will not rerun if the expected output file designated by the return statement in it's `output` method already exists.
