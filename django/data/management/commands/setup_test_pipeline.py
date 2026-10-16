"""
Create the test_pipeline schema as a structural clone of the pipeline schema.

Usage:
    python manage.py setup_test_pipeline

Drops and recreates test_pipeline with the same table structure (columns,
indexes, constraints, sequences) as the pipeline schema, but with no data.
Intended for use in verification testing of pipeline output.
"""

import os
import re
import subprocess

from django.core.management.base import BaseCommand
from django.db import connections


class Command(BaseCommand):
    help = 'Create test_pipeline schema as a structural clone of pipeline'

    def handle(self, *args, **options):
        s = connections['pipeline'].settings_dict

        env = os.environ.copy()
        env['PGPASSWORD'] = s['PASSWORD']

        pg_args = ['-h', s['HOST'], '-U', s['USER'], s['NAME']]

        self.stdout.write('Dumping pipeline schema structure ...')
        dump = subprocess.run(
            ['pg_dump', '--schema=pipeline', '--schema-only'] + pg_args,
            capture_output=True, text=True, env=env,
        )
        if dump.returncode != 0:
            self.stderr.write(f'pg_dump failed:\n{dump.stderr}')
            raise SystemExit(1)

        ddl = dump.stdout
        ddl = re.sub(r'\bpipeline\.', 'test_pipeline.', ddl)
        ddl = re.sub(r'SET search_path = pipeline\b', 'SET search_path = test_pipeline', ddl)

        preamble = (
            'DROP SCHEMA IF EXISTS test_pipeline CASCADE;\n'
            'CREATE SCHEMA test_pipeline;\n'
        )

        self.stdout.write('Creating test_pipeline schema ...')
        psql = subprocess.run(
            ['psql'] + pg_args,
            input=preamble + ddl, text=True, env=env, capture_output=True,
        )
        if psql.returncode != 0:
            self.stderr.write(f'psql failed:\n{psql.stderr}')
            raise SystemExit(1)

        self.stdout.write(self.style.SUCCESS('test_pipeline schema created successfully'))
