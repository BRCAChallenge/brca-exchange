Dockerization of the `slurm.annotate` script from the VICTOR pipeline, which annotates a VCF file with BayesDel scores (see http://fengbj-laboratory.org/victor-2018-12-06/Manual.html)


## Docker Commands

```
docker build . -t brcachallenge/victor:0.1
docker push brcachallenge/victor:0.1
```

## Data Management

The VICTOR pipeline has fairly large data dependencies. In order to allow for reproducible pipeline runs the complete data is redownloaded for every brcaexchange pipeline release. Since large parts of the data are not expected to change much over time, harddisk space can be significantly reduced by hardlinking files which are identical across different releases. This can easily be done with tools like `jdupes` (https://github.com/jbruchon/jdupes).
