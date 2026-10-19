#!/usr/bin/env python
"""
Deduplicate exLOVD VCF records that share the same VRS digest.

When multiple records map to the same variant (same ALT-allele VRS digest),
keep the one whose key_observational_reference cites the most recent year.
Ties are broken by preferring the record with more non-missing INFO fields.

Usage:
    python deduplicate_exlovd_vcf.py --input-vcf IN.vcf.gz --output-vcf OUT.vcf.gz
"""

import argparse
import re
import sys
from collections import defaultdict

import pysam


def _most_recent_year(comments):
    """Return the largest 4-digit year found in the comments value, or 0."""
    if not comments:
        return 0
    text = ' '.join(str(c) for c in comments) if isinstance(comments, (tuple, list)) else str(comments)
    years = [int(y) for y in re.findall(r'\b(?:19|20)\d{2}\b', text)]
    return max(years) if years else 0


def _non_missing_count(rec):
    """Count INFO fields whose value is not '-' (a stand-in for missing)."""
    count = 0
    for key in rec.info:
        val = rec.info[key]
        vals = val if isinstance(val, (tuple, list)) else (val,)
        if any(str(v) != '-' for v in vals):
            count += 1
    return count


def _alt_digest(rec):
    """Return the ALT-allele VRS digest (index 1 of VRS_Allele_IDs), or None."""
    ids = rec.info.get('VRS_Allele_IDs')
    if ids and len(ids) >= 2:
        return ids[1]
    return None


def deduplicate(input_vcf, output_vcf):
    reader = pysam.VariantFile(input_vcf, 'r')
    all_records = [(r.copy(), _alt_digest(r)) for r in reader]
    reader.close()

    by_digest = defaultdict(list)
    for idx, (_, digest) in enumerate(all_records):
        if digest:
            by_digest[digest].append(idx)

    keep = set()
    n_dup_groups = 0

    for digest, indices in by_digest.items():
        if len(indices) == 1:
            keep.add(indices[0])
        else:
            n_dup_groups += 1
            best = max(
                indices,
                key=lambda i: (
                    _most_recent_year(all_records[i][0].info.get('key_observational_reference')),
                    _non_missing_count(all_records[i][0]),
                ),
            )
            keep.add(best)

    for idx, (_, digest) in enumerate(all_records):
        if digest is None:
            keep.add(idx)

    with pysam.VariantFile(input_vcf, 'r') as hdr_src:
        with pysam.VariantFile(output_vcf, 'wz', header=hdr_src.header) as writer:
            for idx, (rec, _) in enumerate(all_records):
                if idx in keep:
                    writer.write(rec)

    n_in = len(all_records)
    n_out = len(keep)
    print(f'exLOVD deduplication: {n_in} records in, {n_out} out, '
          f'{n_dup_groups} duplicate groups resolved, {n_in - n_out} records removed',
          file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input-vcf',  required=True, help='VRS-annotated exLOVD VCF (.vcf.gz)')
    parser.add_argument('--output-vcf', required=True, help='Deduplicated output VCF (.vcf.gz)')
    args = parser.parse_args()
    deduplicate(args.input_vcf, args.output_vcf)


if __name__ == '__main__':
    main()
