# Generated manually on 2026-10-02

from django.db import migrations

# end_pos was always declared NOT NULL by the model/migration 0001, but the
# live column allowed NULL and every row was unbackfilled. Compute it from
# the already-correct pos/ref using the same formula the pipeline uses
# elsewhere (query_clingen_allele_registry.py, load_vcf.py):
# end_pos = pos + len(ref) - 1, preserving the '-' sentinel for rows with no
# valid coordinate on this assembly.
BACKFILL_AND_CONSTRAIN_SQL = """
UPDATE variant_genomic_coordinates
SET end_pos = CASE
    WHEN pos = '-' THEN '-'
    ELSE (pos::integer + length(ref) - 1)::text
END
WHERE end_pos IS NULL;

ALTER TABLE variant_genomic_coordinates ALTER COLUMN end_pos SET NOT NULL;
"""

REVERSE_SQL = """
ALTER TABLE variant_genomic_coordinates ALTER COLUMN end_pos DROP NOT NULL;
"""


class Migration(migrations.Migration):

    dependencies = [
        ('data', '0005_variant_clinvar_aggregate'),
    ]

    operations = [
        migrations.RunSQL(
            sql=BACKFILL_AND_CONSTRAIN_SQL,
            reverse_sql=REVERSE_SQL,
            state_operations=[],
        ),
    ]
