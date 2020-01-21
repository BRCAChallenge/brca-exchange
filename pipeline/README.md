# How to Run "The Pipeline"

## Introduction

In order to faciliate working with the pipeline a makefile is included in the code base under `pipeline/Makefile`. It is parametrized using a configuration file which by default is assumed to be in `pipeline/brca_pipeline_cfg.mk`.

### Requirements

In order to successfully run the pipeline using this setup your environment needs to be equipped with the following:

* docker
* GNU make >=v3.82
* python with jinja2-cli installed (`pip install jinja2-cli`)

#### Mac OS X

Additional steps have to be taken on Mac OS X:

* brew install make --with-default-names
* brew install coreutils
* setup bash version >= 4.0 (https://coderwall.com/p/dmuxma/upgrade-bash-on-your-mac-os)

## Creating a New Data Release

New data releases should ideally be generated on a dedicated pipeline machine. Although it could be run on any other machine in principle, some sources (e.g. LOVD) are only available from there.

### Create a Data Release
To create a new data release the entry point is the `pipeline/pipeline_running/generate_release.sh` script. It needs to be invoked with appropriate arguments, i.e.

 * root working directory
 * path to luigi credentials file (see section below)
 * directory where previous release archives are stored

For the pipeline machine, we for example get:

```
/home/brca/brca_upstream/pipeline/pipeline_running/generate_release.sh /data/monthly_releases /data/luigi_pipeline_credentials.cfg /data/previous_releases
```

This script clones the BRCA Exchange repo into a directory in `WORKING_DIR/data_release_yyyy-MM-dd` referring to the current date and checks out the latest commit on master. It then generates a configuration file `brca_pipeline_cfg.mk` where paths and other settings are set up.

Finally, the following steps are done via the Makefile:
 * downloads resources files
 * builds a docker image
 * kicks off the pipeline in the docker image just created

Should anything go wrong, the pipeline can be easily restarted by issuing `make run-pipeline` in the `pipeline` directory of the code base of the corresponding release (that's where both the `Makefile` and the configuration in `brca_pipeline_cfg.mk` is stored).

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
│   ├── brca_out                <-- pipeline working directory
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
