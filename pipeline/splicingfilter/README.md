Splicing Data Filter
---

For variants listed in a blacklist, replaces in-silico pathogenicity values computed in the 'splicing' step with '-'.

The list of blacklisted variants is contained in `blacklisted_vars.txt`. The variant identifiers consist of the
variant's reference cDNA sequence followed by a ':' and then its hg38 genomic coordinate.

The filtering script, `filterBlacklistedVars.py`, takes the following arguments:
- `--blacklisted_vars FILENAME`: File of HGVS cDNA coordinates of variants which should return blank priors data.
- `--output FILENAME`: File to write filtered variants, defaults to stdout.
- `--use-csv BOOLEAN`: Use csv module (instead of pandas) to transform input, which is slower but more accurate.

When pandas reads in a DataFrame from a CSV via `pd.read_csv()` it parse the float fields, introducing some imprecision
due to rounding errors. The errors are small, but they're sufficient to change the checksum and equality checks in the
tests, which is why the csv-module implementation is used as the default (despite being slightly slower), and why the
pandas test is marked as expected to fail.