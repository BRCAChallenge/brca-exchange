

from django.db import migrations


class Migration(migrations.Migration):
    dependencies = [
        ('data', '0035_add_functional_analysis_fields'),
    ]

    operations = [
        migrations.RunSQL("""
        CREATE INDEX report_submission_id_lovd ON report ("Submission_ID_LOVD");
    """)

    ]
