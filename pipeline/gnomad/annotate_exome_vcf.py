#!/usr/bin/env python
"""One-off script: strip pre-existing VRS fields from exome VCF and run vrs-annotate."""

import os
import re
import subprocess
import sys
import tempfile

EXOME_IN  = "/data/new_schema/resources/gnomADv4.1.exome.hg38.vcf"
EXOME_VCF = "/data/new_schema/data_out/vcf/gnomADv4.1.exome.hg38.vcf.gz"
EXOME_PKL = "/data/new_schema/data_out/vcf/gnomADv4.1.exome.hg38.dicts.pkl"
SEQREPO   = "/data/seqrepo/latest"

_GRCh38_CONTIGS = {
    'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
    'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636,
    'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
    'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
    'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415,
}


def normalize_vcf_header(input_vcf, output_vcf):
    header, data = [], []
    print(f"Reading {input_vcf} ...", flush=True)
    with open(input_vcf) as fh:
        for line in fh:
            if line.startswith('#'):
                header.append(line)
            else:
                data.append(line)

    declared_contigs = set()
    for line in header:
        m = re.search(r'^##contig=<ID=([^,>]+)', line)
        if m:
            declared_contigs.add(m.group(1))
    referenced_contigs = {row.split('\t')[0] for row in data if row.strip()}
    missing_contigs = [c for c in (referenced_contigs - declared_contigs)
                       if c in _GRCh38_CONTIGS]

    declared_info = set()
    for line in header:
        m = re.search(r'^##INFO=<ID=([^,>]+)', line)
        if m:
            declared_info.add(m.group(1))
    referenced_info = set()
    for row in data:
        if not row.strip():
            continue
        cols = row.split('\t')
        if len(cols) > 7 and cols[7] not in ('.', ''):
            for field in cols[7].split(';'):
                key = field.split('=')[0]
                if key:
                    referenced_info.add(key)
    missing_info = sorted(referenced_info - declared_info)

    print(f"  missing_contigs={missing_contigs}, missing_info={missing_info}", flush=True)

    if not missing_contigs and not missing_info:
        with open(output_vcf, 'w') as fh:
            fh.writelines(header)
            fh.writelines(data)
        return

    chrom_idx = next(i for i, l in enumerate(header) if l.startswith('#CHROM'))
    inject = []
    if missing_contigs:
        inject += [f'##contig=<ID={c},length={_GRCh38_CONTIGS[c]},assembly=GRCh38>\n'
                   for c in sorted(missing_contigs)]
    if missing_info:
        inject += [f'##INFO=<ID={k},Number=.,Type=String,Description="{k}">\n'
                   for k in missing_info]
    header = header[:chrom_idx] + inject + header[chrom_idx:]
    with open(output_vcf, 'w') as fh:
        fh.writelines(header)
        fh.writelines(data)
    print(f"  Wrote normalized VCF to {output_vcf}", flush=True)


def run(cmd):
    print(f"Running: {' '.join(cmd)}", flush=True)
    result = subprocess.run(cmd, check=True)
    return result


def main():
    with tempfile.NamedTemporaryFile(suffix='_stripped.vcf', delete=False) as t1:
        stripped = t1.name
    with tempfile.NamedTemporaryFile(suffix='_normalized.vcf', delete=False) as t2:
        normalized = t2.name

    try:
        print("Step 1: Strip pre-existing VRS INFO fields ...", flush=True)
        run([
            'bcftools', 'annotate', '-x',
            'INFO/VRS_Allele_IDs,INFO/VRS_Starts,INFO/VRS_Ends,INFO/VRS_States',
            '-o', stripped, EXOME_IN
        ])
        print(f"  Stripped VCF: {os.path.getsize(stripped) // (1024*1024)} MB", flush=True)

        print("Step 2: Normalize VCF header ...", flush=True)
        normalize_vcf_header(stripped, normalized)
        print(f"  Normalized VCF: {os.path.getsize(normalized) // (1024*1024)} MB", flush=True)

        print("Step 3: VRS-annotate ...", flush=True)
        run([
            'vrs-annotate', 'vcf', normalized,
            '--vcf-out', EXOME_VCF,
            '--pkl-out', EXOME_PKL,
            '--dataproxy-uri', f'seqrepo+file://{SEQREPO}',
        ])
        print("Done.", flush=True)
        print(f"  {EXOME_VCF}: {os.path.getsize(EXOME_VCF) // (1024*1024)} MB", flush=True)
        print(f"  {EXOME_PKL}: {os.path.getsize(EXOME_PKL) // (1024*1024)} MB", flush=True)

    finally:
        for f in (stripped, normalized):
            try:
                os.unlink(f)
            except FileNotFoundError:
                pass


if __name__ == '__main__':
    main()
