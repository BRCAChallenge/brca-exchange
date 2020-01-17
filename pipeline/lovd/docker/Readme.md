ExtractDataFromLatestEXLOVD was dockerized in order not to have to port the
`extract_data.py` script to python 3.

## Docker Commands

```
docker build . -t brcachallenge/exlovd-extraction:0.1
docker push brcachallenge/exlovd-extraction:0.1
```
