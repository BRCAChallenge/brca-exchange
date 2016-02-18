from __future__ import unicode_literals

from django.db import migrations


class Migration(migrations.Migration):
    dependencies = [
        ('data', '0003_data'),
    ]

    operations = [
        migrations.RunSQL("""
        CREATE EXTENSION pg_trgm;

        CREATE TABLE words AS SELECT word FROM ts_stat('SELECT fts_document FROM variant');

        CREATE INDEX words_idx ON words USING gin(word gin_trgm_ops);
    """)

    ]
