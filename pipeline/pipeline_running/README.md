# How to Run "The Pipeline"

## Introduction

In order to faciliate working with the pipeline a makefile is included in the code base under `pipeline/Makefile`. It is parametrized using a configuration file which by default is assumed to be in `pipeline/brca_pipeline_cfg.mk`.

## Creating a New Data Release

New data releases should ideally be generated on a dedicated pipeline machine. Although it could be run on any other machine in principle, some sources (e.g. LOVD) are only available from the pipeline machine.

### Credentials:

Credentials can be passed into the container by mounting an appropriate file. Currently, the file should contain the following:

```
[RunAll]
# BIC credentials
u=bicusername
p=bicpassword

# synapse credentials
synapse_username=your_username
synapse_password=some_password
synapse_enigma_file_id=syn8465585

```


### Create the data
To create a new data release entry point is the `pipeline/pipeline_running/generate_release.sh` scripts. It needs to be started with appropriate arguments, i.e.

 * root work directory
 * path to luigi credentials file (see above)
 * directory where previous release archives are stored

For the pipeline machine, we get:

```
./generate_release.sh /home/pipeline/monthly_releases /home/pipeline/luigi_pipeline_credentials.cfg /home/pipeline/previous_releases
```

This script clones the BRCA Exchange repo into an appropriate directory and checks out the latest commit on master. It then generates a confguration file `brca_pipeline_cfg.mk` where paths and other settings are set up.

Finally, the following steps are done via the Makefile: 
 * generates appropriate directories, every release happens in a separate directory. A release is named `data_release_yyyy-MM-dd` referring to the current date.
 * checks out the git repository code on the commit hash specified above
 * downloads resources files
 * builds a docker image
 * creates a file with commands to be run after the release data has been validated
 * kicks off the pipeline in the docker image just created

Should anything go wrong, the pipeline can easily be restarted by issuing `make run-pipeline` in the `pipeline` directory of the code base.

### Postprocessing

After the data in the tar release file has been sanity checked (and the release notes updated), some post processing steps need to be done.

Steps include:
 * updating the release notes and regenerating the release archive with the release notes
 * tagging the commit in the main git repository
 * pushing the docker image to dockerhub
 * copying the release tar to `previous_releases` folder.

This can be done in one breeze by running `make post-release-cmds`.

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

TODO: coming soon
