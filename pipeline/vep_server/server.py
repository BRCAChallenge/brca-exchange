"""
Minimal VEP REST server.

Endpoints:
  GET /vep/hgvs/<hgvs_notation>  — run VEP on one HGVS string, return JSON array
  GET /health                     — liveness check
"""

import json
import os
import subprocess
import tempfile
from http.server import BaseHTTPRequestHandler, HTTPServer
from urllib.parse import unquote

VEP_CACHE_DIR = os.environ.get('VEP_CACHE_DIR', '/vep_cache')
ASSEMBLY      = os.environ.get('VEP_ASSEMBLY',   'GRCh38')
PORT          = int(os.environ.get('VEP_PORT',   '8080'))


def run_vep(hgvs_notation):
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as fh:
        fh.write(hgvs_notation + '\n')
        tmp = fh.name
    try:
        result = subprocess.run(
            [
                'vep',
                '--input_file',  tmp,
                '--format',      'hgvs',
                '--output_file', 'STDOUT',
                '--json',
                '--no_stats',
                '--cache',
                '--dir_cache',   VEP_CACHE_DIR,
                '--offline',
                '--assembly',    ASSEMBLY,
                '--everything',
            ],
            capture_output=True, text=True, timeout=120,
        )
        if result.returncode != 0:
            raise RuntimeError(result.stderr.strip())
        records = [
            json.loads(line)
            for line in result.stdout.splitlines()
            if line.strip()
        ]
        return records
    finally:
        os.unlink(tmp)


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


if __name__ == '__main__':
    print(f'VEP server starting on port {PORT} (assembly={ASSEMBLY}, cache={VEP_CACHE_DIR})',
          flush=True)
    HTTPServer(('0.0.0.0', PORT), Handler).serve_forever()
