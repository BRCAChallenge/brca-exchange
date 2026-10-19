"""
Tests for compare_table.py: validate variant_exlovd in pipeline vs test_pipeline.

Connects to the live database (storage.pg) and asserts structural invariants:
  - all rows from test_pipeline appear in pipeline (deleted = 0)
  - no rows in pipeline are absent from test_pipeline (added = 0)
  - shared row count equals the expected total

Run with:
    pytest diff/test_compare_table.py
"""

import tempfile
import pytest
import psycopg2
import psycopg2.extras

from compare_table import get_pk_columns, get_columns

DB = dict(host='localhost', dbname='storage.pg', user='postgres', password='postgres')
OLD = 'test_pipeline'
NEW = 'pipeline'
TABLE = 'variant_exlovd'
EXPECTED_ROW_COUNT = 1628


@pytest.fixture(scope='module')
def conn():
    c = psycopg2.connect(**DB)
    c.autocommit = True
    yield c
    c.close()


@pytest.fixture(scope='module')
def cur(conn):
    with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as c:
        yield c


def test_pk_is_vrs_digest(cur):
    pk = get_pk_columns(cur, OLD, TABLE)
    assert pk == ['VRS_Digest']


def test_columns_match_across_schemas(cur):
    old_cols = get_columns(cur, OLD, TABLE)
    new_cols = get_columns(cur, NEW, TABLE)
    assert old_cols == new_cols, (
        f"Column mismatch between {OLD}.{TABLE} and {NEW}.{TABLE}: "
        f"old={old_cols}, new={new_cols}"
    )


def test_no_rows_deleted(cur):
    """Every row in test_pipeline must appear in pipeline (same VRS_Digest)."""
    cur.execute(f"""
        SELECT COUNT(*) FROM "{OLD}"."{TABLE}" o
        LEFT JOIN "{NEW}"."{TABLE}" n ON n."VRS_Digest" = o."VRS_Digest"
        WHERE n."VRS_Digest" IS NULL
    """)
    deleted = cur.fetchone()['count']
    assert deleted == 0, f"{deleted} rows in {OLD}.{TABLE} are absent from {NEW}.{TABLE}"


def test_no_rows_added(cur):
    """Every row in pipeline must appear in test_pipeline (same VRS_Digest)."""
    cur.execute(f"""
        SELECT COUNT(*) FROM "{NEW}"."{TABLE}" n
        LEFT JOIN "{OLD}"."{TABLE}" o ON o."VRS_Digest" = n."VRS_Digest"
        WHERE o."VRS_Digest" IS NULL
    """)
    added = cur.fetchone()['count']
    assert added == 0, f"{added} rows in {NEW}.{TABLE} are absent from {OLD}.{TABLE}"


def test_shared_row_count(cur):
    cur.execute(f"""
        SELECT COUNT(*) FROM "{OLD}"."{TABLE}" o
        JOIN "{NEW}"."{TABLE}" n ON n."VRS_Digest" = o."VRS_Digest"
    """)
    shared = cur.fetchone()['count']
    assert shared == EXPECTED_ROW_COUNT, (
        f"Expected {EXPECTED_ROW_COUNT} shared rows, got {shared}"
    )


def test_diff_output_is_written():
    """End-to-end: compare_table.py runs without error and writes output."""
    import subprocess
    import sys
    import os

    script = os.path.join(os.path.dirname(__file__), 'compare_table.py')
    with tempfile.NamedTemporaryFile(suffix='.txt', delete=False) as fh:
        out_path = fh.name

    result = subprocess.run(
        [sys.executable, script,
         '--table', TABLE,
         '--old', OLD,
         '--new', NEW,
         '--output', out_path],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"compare_table.py failed:\n{result.stderr}"
    assert 'shared:' in result.stdout
    assert 'deleted:' in result.stdout
    assert 'added:' in result.stdout
    assert os.path.exists(out_path) and os.path.getsize(out_path) > 0
