# Luigi (WIP)

Contains code used to stitch together processes related to the pipeline.

## Instructions

To run, you must first install Luigi: `pip install luigi`.

Also, make sure that all dependencies of the pipeline directory are properly installed (see https://github.com/BD2KGenomics/brca-exchange/tree/luigi/pipeline).

Create environment variables for paths to an output directory, BRCA resources directory, and a file parent directory. The file parent directory will hold directories for each of the separate data sources after running the script.

Username and password are required to download files from BIC. They can be found in `/hive/groups/cgl/brca/phase1/data/bic/account.txt` at UCSC.

Synapse username and password are also required to download files for Enigma. They can be created at [synapse.org](http://synapse.org). You will also need access to the Enigma directory with Synapse ID `syn7188267` and the id for the preprocessed Enigma combined output file.

The Luigi script takes several arguments to run:

* `--module`: names the module containing the luigi script to run followed by the name of the luigi task to run (to run all tasks, run the master task)
* `--u`: username for access to BIC data
* `--p`: password for access to BIC data
* `--synapse-username`: username from synapse.org to download Enigma data
* `--synapse-password`: password from synapse.org to download Enigma data
* `--synapse-enigma-file-id`: id of the preprocessed Enigma combined output file in synapse
* `--output-dir`: the directory to store output from all sources in the luigi script
* `--resources-dir`: the directory containing BRCA Resources obtained by following instructions in the pipeline readme
* `--file-parent-dir`: the directory to store output from individual sources in the luigi script
* `--previous-release` (optional): the previous data release used to compare with the newest release to determine variant changes between releases. If this argument is not provided, changes to variants between release versions will not be determined or appended to the output file.

To run: `python -m luigi --module CompileVCFFiles RunAll --u {username} --p {password} --synapse-username {username from synapse.org} --synapse-password {password from synapse.org} --synapse-enigma-file-id {id for combined enigma output file from synapse} --output-dir $OUTPUT_DIR --resources-dir $BRCA_RESOURCES --file-parent-dir $PARENT_DIR --previous-release $PREVIOUS_RELEASE --local-scheduler`

You can replace `RunAll` with individual tasks, or comment out required tasks in the `RunAll` task to control which tasks are run. A task will not rerun if the expected output file designated by the return statement in it's `output` method already exists.
