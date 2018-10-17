Splicing Data Filter
---

For variants listed in a blacklist, replaces in-silico pathogenicity values computed in the 'splicing' step with '-'.

The list of blacklisted variants is contained in `blacklisted_vars.txt`. The variant identifiers consist of the
variant's reference cDNA sequence followed by a ':' and then its hg38 genomic coordinate.

The filtering script, `filterBlacklistedVars.py`, takes the following arguments:
- `--blacklisted_vars FILENAME`: File of HGVS cDNA coordinates of variants which should return blank priors data.
- `--output FILENAME`: File to write filtered variants, defaults to stdout.

When pandas reads in a DataFrame from a CSV via `pd.read_csv()` it parses the float fields, introducing some imprecision
due to rounding errors. The test accommodates this imprecision by checking if float-coercible values are almost equal to
each other, using absolute equality for the other values.