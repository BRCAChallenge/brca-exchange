# -*- coding: utf-8 -*-


from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('data', '0014_add_lovd_fields'),
    ]

    operations = [
        migrations.AddField(
            model_name='report',
            name='AA_Allele_Frequency_ESP',
            field=models.TextField(default='-'),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='report',
            name='Allele_Frequency_ESP',
            field=models.TextField(default='-'),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='report',
            name='EA_Allele_Frequency_ESP',
            field=models.TextField(default='-'),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='variant',
            name='AA_Allele_Frequency_ESP',
            field=models.TextField(default='-'),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='variant',
            name='Allele_Frequency_ESP',
            field=models.TextField(default='-'),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='variant',
            name='EA_Allele_Frequency_ESP',
            field=models.TextField(default='-'),
            preserve_default=False,
        ),
        migrations.RunSQL(
            """
            DROP MATERIALIZED VIEW IF EXISTS currentvariant;
            CREATE MATERIALIZED VIEW currentvariant AS (
                SELECT * FROM "variant" WHERE (
                    "id" IN ( SELECT DISTINCT ON ("Genomic_Coordinate_hg38") "id" FROM "variant" ORDER BY "Genomic_Coordinate_hg38" ASC, "Data_Release_id" DESC )
                )
            );
            """
        ),
    ]
