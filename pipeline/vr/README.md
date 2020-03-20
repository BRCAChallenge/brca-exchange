appendVRIds was dockerized in order to run using py3 (https://github.com/ga4gh/vr-python requires python 3) while the existing pipeline is still running on py2

## Docker Commands

```
docker build . -t brcachallenge/append-vr-ids:0.1
docker push brcachallenge/append-vr-ids:0.1
```