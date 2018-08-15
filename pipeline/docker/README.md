# BRCA Exchange Docker Pipeline

## Build Image

From the `pipeline/docker/` directory, run the script `build_docker_image.sh` to build a docker image on your system.

## Running Container

### Preparations

Various files and directories need to be prepared and linked into the container in which the pipeline is run:

 * Resources (e.g. genome annotations) can be downloaded into a directory using `pipeline/download_resource_files.sh`
 * Output directory to write the results to
 * code base of brca-exchange (see below)
 * files containing credentials
 * previous release files
 * release notes for the currently generated release

A version of the code base is already included in the docker image. To use a different version, e.g. during development or to reproduce a result from a specific version, the code base can be easily replaced by mounting some appropriate directory on the host system.

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

### Running the Pipeline
Below an example invocation of `docker run`. In the first line we make sure, that the container is run as the same user as the one invoking this command instead of `root`. This avoids issues further downstream, where the output file would be owned by `root`.

Make sure the directories on the host system already exist (in particular the corresponding directory of `/files/data`), as otherwise issues may arise with access rights.

In the following lines, paths on the host are mapped to paths in the container. You would need to adapt the path before the `:` accordingly.

Notes:

* The line concerning the code base can be omitted. In this case, the version of the pipeline already contained within the image is run.
* The line concerning the DATA_DATE can be omitted. If provided, it explicitly sets the date of the pipeline run rather than using the current date (e.g. this is useful when rerunning a pipeline that aborted half way through).
* The line concerning UTA_DB_URL can be used to run a local instance of a uta database. This is helpful if a remote instance is having issues.

```
docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) \
       -e "DATA_DATE=${DATA_DATE}" \
       --network host \
       -e "UTA_DB_URL=postgresql://anonymous@0.0.0.0:50827/uta/uta_20170629" \
       -v  path_to_resource_files:/files/resources \
       -v  path_to_output_directory:/files/data \
       -v  optional_path_to_code_base:/opt/brca-exchange \
       -v  path_to_pipeline_credentials.cfg:/opt/luigi_pipeline_credentials.cfg \
       -v  path_to_previous_release_tar:/files/previous_release.tar.gz \
       -v  path_to_release_notes.txt:/files/release_notes.txt \
       -v /tmp:/.synapseCache \
       -it \
       brcachallenge/brca-exchange-pipeline
```
