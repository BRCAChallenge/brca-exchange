"""
Minimal VEP REST server.

VEP forbids HGVS format in offline mode, so genomic HGVS strings are
converted to VCF format internally before being passed to VEP.

Endpoints:
  GET /vep/hgvs/<hgvs_notation>  — run VEP on one genomic HGVS string, return JSON array
  GET /health                     — liveness check
"""

import json
import os
import re
import subprocess
import tempfile
from http.server import BaseHTTPRequestHandler, HTTPServer
from urllib.parse import unquote

VEP_CACHE_DIR = os.environ.get('VEP_CACHE_DIR', '/vep_cache')
VEP_FASTA     = os.environ.get('VEP_FASTA',     '/references/hg38.fa')
ASSEMBLY      = os.environ.get('VEP_ASSEMBLY',   'GRCh38')
PORT          = int(os.environ.get('VEP_PORT',   '8080'))

# Map RefSeq contig accessions to UCSC-style chromosome names (chr-prefixed)
CONTIG_TO_CHROM = {
    'NC_000001.11': 'chr1',  'NC_000002.12': 'chr2',  'NC_000003.12': 'chr3',
    'NC_000004.12': 'chr4',  'NC_000005.10': 'chr5',  'NC_000006.12': 'chr6',
    'NC_000007.14': 'chr7',  'NC_000008.11': 'chr8',  'NC_000009.12': 'chr9',
    'NC_000010.11': 'chr10', 'NC_000011.10': 'chr11', 'NC_000012.12': 'chr12',
    'NC_000013.11': 'chr13', 'NC_000014.9':  'chr14', 'NC_000015.10': 'chr15',
    'NC_000016.10': 'chr16', 'NC_000017.11': 'chr17', 'NC_000018.10': 'chr18',
    'NC_000019.10': 'chr19', 'NC_000020.11': 'chr20', 'NC_000021.9':  'chr21',
    'NC_000022.11': 'chr22', 'NC_000023.11': 'chrX',  'NC_000024.10': 'chrY',
    'NC_012920.1':  'chrM',
}


def _fasta_fetch(chrom, start, end):
    """Return the reference sequence for chrom:start-end (1-based, inclusive)."""
    result = subprocess.run(
        ['samtools', 'faidx', VEP_FASTA, f'{chrom}:{start}-{end}'],
        capture_output=True, text=True,
    )
    lines = result.stdout.splitlines()
    return ''.join(lines[1:]).upper()


def genomic_hgvs_to_vcf(hgvs_str):
    """
    Convert a genomic HGVS string to a VCF record tuple (chrom, pos, ref, alt).

    Handles: SNV, MNV, del, ins, delins, dup, identity (=).
    Raises ValueError for unrecognised patterns.
    """
    m = re.match(r'(NC_\d+\.\d+):g\.(.+)', hgvs_str)
    if not m:
        raise ValueError(f'Unrecognised HGVS: {hgvs_str}')
    contig, change = m.group(1), m.group(2)
    chrom = CONTIG_TO_CHROM.get(contig)
    if not chrom:
        raise ValueError(f'Unknown contig: {contig}')

    # SNV / MNV:  43094304T>C  or  43094304_43094306ACG>TTT
    sub = re.match(r'(\d+)(?:_\d+)?([ACGT]+)>([ACGT]+)$', change)
    if sub:
        return chrom, int(sub.group(1)), sub.group(2), sub.group(3)

    # Identity:  32355250=
    ident = re.match(r'(\d+)=$', change)
    if ident:
        pos = int(ident.group(1))
        ref = _fasta_fetch(chrom, pos, pos)
        return chrom, pos, ref, ref

    # Deletion:  32315648_32315650del  or  32315648del
    dele = re.match(r'(\d+)(?:_(\d+))?del$', change)
    if dele:
        start = int(dele.group(1))
        end   = int(dele.group(2)) if dele.group(2) else start
        anchor_ref = _fasta_fetch(chrom, start - 1, end)   # includes anchor base
        anchor_alt = anchor_ref[0]
        return chrom, start - 1, anchor_ref, anchor_alt

    # Insertion:  32315648_32315649insACGT
    ins = re.match(r'(\d+)_(\d+)ins([ACGT]+)$', change)
    if ins:
        pos     = int(ins.group(1))
        inserted = ins.group(3)
        anchor  = _fasta_fetch(chrom, pos, pos)
        return chrom, pos, anchor, anchor + inserted

    # Delins:  32315648_32315650delinsACGT  or  32315648delinsACGT
    delins = re.match(r'(\d+)(?:_(\d+))?delins([ACGT]+)$', change)
    if delins:
        start   = int(delins.group(1))
        end     = int(delins.group(2)) if delins.group(2) else start
        ref     = _fasta_fetch(chrom, start, end)
        alt     = delins.group(3)
        return chrom, start, ref, alt

    # Duplication:  32343971_32349243dup
    dup = re.match(r'(\d+)_(\d+)dup$', change)
    if dup:
        start = int(dup.group(1))
        end   = int(dup.group(2))
        seq   = _fasta_fetch(chrom, start, end)
        anchor = _fasta_fetch(chrom, start - 1, start - 1)
        return chrom, start - 1, anchor, anchor + seq

    raise ValueError(f'Unhandled HGVS change: {change}')


VEP_ARGS = [
    'vep',
    '--format',      'vcf',
    '--output_file', 'STDOUT',
    '--json',
    '--no_stats',
    '--cache',
    '--dir_cache',   VEP_CACHE_DIR,
    '--fasta',       VEP_FASTA,
    '--offline',
    '--assembly',    ASSEMBLY,
    '--everything',
]


def _run_vep_on_vcf(vcf_lines):
    """Write vcf_lines to a temp file, run VEP, return parsed JSON records."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as fh:
        fh.write('##fileformat=VCFv4.1\n')
        fh.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        fh.writelines(vcf_lines)
        tmp = fh.name
    try:
        result = subprocess.run(
            VEP_ARGS + ['--input_file', tmp],
            capture_output=True, text=True, timeout=600,
        )
        if result.returncode != 0:
            raise RuntimeError(result.stderr.strip())
        return [json.loads(l) for l in result.stdout.splitlines() if l.strip()]
    finally:
        os.unlink(tmp)


def run_vep(hgvs_notation):
    """Run VEP on a single HGVS string; return list of result records."""
    chrom, pos, ref, alt = genomic_hgvs_to_vcf(hgvs_notation)
    return _run_vep_on_vcf([f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n'])


def run_vep_batch(hgvs_list):
    """Run VEP on a list of HGVS strings in a single VEP call.

    Returns a dict mapping each input HGVS string to its list of VEP records.
    Entries that fail HGVS→VCF conversion are returned as {"error": "..."}.
    """
    vcf_lines = []
    id_to_hgvs = {}   # VCF ID field → original HGVS (for result matching)
    errors = {}

    for i, hgvs in enumerate(hgvs_list):
        vid = f'v{i}'
        try:
            chrom, pos, ref, alt = genomic_hgvs_to_vcf(hgvs)
            vcf_lines.append(f'{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\t.\t.\n')
            id_to_hgvs[vid] = hgvs
        except Exception as e:
            errors[hgvs] = {'error': str(e)}

    results = {hgvs: [] for hgvs in hgvs_list}
    results.update(errors)

    if vcf_lines:
        for record in _run_vep_on_vcf(vcf_lines):
            vid = record.get('id', '')
            hgvs = id_to_hgvs.get(vid)
            if hgvs:
                results[hgvs].append(record)

    return results


class Handler(BaseHTTPRequestHandler):
    def log_message(self, fmt, *args):
        pass  # suppress per-request access log noise

    def _send(self, status, body, content_type='application/json'):
        encoded = body.encode()
        self.send_response(status)
        self.send_header('Content-Type', content_type)
        self.send_header('Content-Length', len(encoded))
        self.end_headers()
        self.wfile.write(encoded)

    def do_GET(self):
        path = unquote(self.path)

        if path == '/health':
            self._send(200, '{"status":"ok"}')
            return

        if path.startswith('/vep/hgvs/'):
            hgvs = path[len('/vep/hgvs/'):]
            if not hgvs:
                self._send(400, '{"error":"missing HGVS notation"}')
                return
            try:
                records = run_vep(hgvs)
                self._send(200, json.dumps(records))
            except Exception as e:
                self._send(500, json.dumps({'error': str(e)}))
            return

        self._send(404, '{"error":"not found"}')

    def do_POST(self):
        path = unquote(self.path)

        if path == '/vep/hgvs/batch':
            length = int(self.headers.get('Content-Length', 0))
            body = self.rfile.read(length)
            try:
                hgvs_list = json.loads(body)
                if not isinstance(hgvs_list, list):
                    raise ValueError('Expected a JSON array')
                results = run_vep_batch(hgvs_list)
                self._send(200, json.dumps(results))
            except Exception as e:
                self._send(400, json.dumps({'error': str(e)}))
            return

        self._send(404, '{"error":"not found"}')


if __name__ == '__main__':
    print(f'VEP server starting on port {PORT} (assembly={ASSEMBLY}, cache={VEP_CACHE_DIR})',
          flush=True)
    HTTPServer(('0.0.0.0', PORT), Handler).serve_forever()
