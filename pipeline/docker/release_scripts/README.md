## Data Release Procedure

### Create the data
To create a new data release the `create_monthly_release.sh` needs to be started with appropriate arguments, i.e. 

 * the git hash tag of the main repo
 * path to previous release tar
 * a release tag (optional)

The release tag is used to tag docker images and the release commits and also file names. By default it is `data_release_yyyy-MM-dd` referring to the current date.

This script does the following: 
 * generates appropriate directories (every release happens in a separate directory)
 * checks out the git repository code on the commit hash specified above
 * downloads resources files
 * builds a docker image
 * creates a file with commands to be run after the release data has been validated
 * kicks off the pipeline in the docker image just created

### Postprocessing

After the data in the tar release file has been sanity checked, some post processing steps need to be done. To facilitate this for every release a file with useful commands is generated in the directory `monthly_releases/scripts/postprocessing_commands`.

Steps include:
 * updating the release notes and regenerating the release archive with the release notes
 * tagging the commit in the main git repository
 * pushing the docker image to dockerhub
 * copying the release tar to `previous_releases` folder.

## Involved Directories

```
monthly_releases
├── data_release_TAG            <-- release working dir
│   ├── brca-exchange           <-- clone of git repository 
│   ├── brca_out                <-- pipeline working directory
│   └── resources               <-- e.g. reference sequences
├── logs                        <-- log data of pipeline runs 
├── release_notes               <-- release notes of the various releases
├── scripts                     <-- scripts and templates to run pipeline
│   └── postprocessing_commands <-- postprocessing commands for various releases
previous_releases               <-- released archives of previous releases
```

## Crontab

To run the pipeline regularly the following crontab is added, where every Friday at 6 in the morning UTC the pipeline is kicked off.

```
MAILTO=PUT_ADDRESS1_HERE,PUT_ADDRESS2_HERE
SHELL=/bin/bash
0 6 * * 5 source /home/pipeline/.profile; /home/pipeline/monthly_releases/scripts/create_monthly_release.sh master /home/pipeline/previous_releases/latest_release.tar.gz data_release_$(date +"\%Y-\%m-\%d")
```
