from __future__ import unicode_literals

from django.db import migrations


class Migration(migrations.Migration):
    dependencies = [
        ('data', '0030_add_custom_lovd_submission_ids'),
    ]

    operations = [
        migrations.RunSQL("""
        CREATE INDEX report_variant_id ON report ("Variant_id");
        CREATE INDEX report_scv_clinvar ON report ("SCV_ClinVar");
        CREATE INDEX report_dbid_lovd ON report ("DBID_LOVD");
        CREATE INDEX report_source ON report ("Source");
    """)

    ]
