# How to Run "The Pipeline"

## Introduction

In order to faciliate working with the pipeline a makefile is included in the code base under `pipeline/Makefile`. It is parametrized using a configuration file which by default is assumed to be in `pipeline/brca_pipeline_cfg.mk`.

### Requirements
* docker
* GNU make >=v3.82
* python with jinja2-cli installed (`pip install jinja2-cli`)

## Data Releases

New data releases should ideally be generated on a dedicated pipeline machine. Although it could be run on any other machine in principle, some sources (e.g. LOVD) are only available from there.

### Create a Data Release

To create a new data release the entry point is the `pipeline/pipeline_running/generate_release.py` script. You must be checked out on a branch that has been pushed to GitHub because the first thing the script does is clone the repo from scratch into the new working directory.

#### Mandatory Arguments
 * `root_dir` : Root directory for the release (working directory will be created here)
 * `credentials_path` : Path to Luigi credentials configuration file

#### Optional Arguments
 * `gene_config_filename` : Gene configuration filename (default: `gene_config_brca_only.txt`)
 * `--previous-release-dir` : Directory containing previous release bundle, for comparison to the new release (default: the new working directory)

#### Usage Example

```
git clone https://github.com/BRCAChallenge/brca-exchange.git
~/brca-exchange/pipeline/pipeline_running/generate_release.py \
	/data/monthly_releases \
	/data/luigi_pipeline_credentials.cfg \
	gene_config_brca_only.txt
```

This script clones the BRCA Exchange repo into a directory `data_release_yyyy-MM-dd` within `root_dir` referring to the current date and checks out the same commit as the repo you invoked the script with. It then generates a configuration file `brca_pipeline_cfg.mk` where paths and other settings are set up.

Then, the following steps are done via the Makefile:
 * downloads resources files
 * builds a docker image
 * kicks off the pipeline in the docker image just created

Should anything go wrong, the pipeline can be easily restarted by issuing `make build-release` from the generated data release dir `data_release_yyyy-MM-dd/code/pipeline` (that's where both the `Makefile` and the configuration in `brca_pipeline_cfg.mk` are stored).

### Postprocessing

After the data in the tar release file has been sanity checked (and the release notes updated), some post processing steps need to be done.

Steps include:
 * updating the release notes and regenerating the release archive with the release notes
 * tagging the commit in the main git repository
 * pushing the docker image to dockerhub
 * copying the release tar to `previous_releases` folder.

This can be done in one breeze by running `make post-release-cmds`.

### Credentials

Early stages of the pipeline need credentials to download data. These can be passed into the container by mounting an appropriate file. Also note, that some data sets are only available via the pipeline machine. However, later stages of the pipeline don't need any and a dummy file could be created.

Currently, such a credential files should contain the following:

```
[PipelineParams]
# BIC credentials
u=bicusername
p=bicpassword

```

### Setup on Pipeline Machine

In directory `/home/pipeline`

```
brca_upstream                   <-- BRCA exchange code base
monthly_releases
├── data_release_TAG            <-- release working dir
│   ├── code                    <-- clone of git repository
│   ├── data_out                <-- pipeline working directory
│   └── resources               <-- e.g. reference sequences
│   └── references              <-- e.g. reference sequences for the splicing pipeline (may be merged in the future)
previous_releases               <-- released archives of previous releases
```

## Developing New Features

A very rough guide on how to use the Makefile target for easier development:

Change to the `pipeline` directory and type the following:

* `make` or `make help` to see what targets are available along with minimal help
* `make init` to set up a configuration file `pipeline/brca_pipeline_cfg.mk` with paths and other settings. It is advisable to edit it according your needs:
* `make setup-dev-env`: runs various targets to set up a dev environment.

Running tasks:
* `make show-luigi-graph`: shows the graph of tasks on the console (use e.g.
`less -R` if you experience issues with colors)
* `make run-interactive`: starts a bash in brca exchange docker container.
* `make run-task [TASK]`: runs a specific luigi task
* `make force-run-task [TASK]`: runs a specific luigi task, deleting its outputs
 first (otherwise luigi doesn't run the task)
* `make clean-files-from [TASK]`: deletes all outputs of the given task along
with all the tasks directly or indirectly depending on it. This is useful to
force regeneration of 'downstream' data if something in TASK has changed.

Running tests:
* `make test`: runs pipeline unit tests in docker container
* `make test-coverage`: runs pipeline unit tests with coverage analysis. The
 HTML reports can be found in the directory `pipline/htmlcov`
