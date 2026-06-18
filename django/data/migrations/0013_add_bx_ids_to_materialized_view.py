# -*- coding: utf-8 -*-


from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ('data', '0012_auto_20170411_0706'),
    ]

    operations = [
        migrations.AddField(
            model_name='CurrentVariant',
            name='BX_ID_1000_Genomes',
            field=models.TextField(default=b''),
        ),
        migrations.AddField(
            model_name='CurrentVariant',
            name='BX_ID_BIC',
            field=models.TextField(default=b''),
        ),
        migrations.AddField(
            model_name='CurrentVariant',
            name='BX_ID_ClinVar',
            field=models.TextField(default=b''),
        ),
        migrations.AddField(
            model_name='CurrentVariant',
            name='BX_ID_ENIGMA',
            field=models.TextField(default=b''),
        ),
        migrations.AddField(
            model_name='CurrentVariant',
            name='BX_ID_ESP',
            field=models.TextField(default=b''),
        ),
        migrations.AddField(
            model_name='CurrentVariant',
            name='BX_ID_ExAC',
            field=models.TextField(default=b''),
        ),
        migrations.AddField(
            model_name='CurrentVariant',
            name='BX_ID_LOVD',
            field=models.TextField(default=b''),
        ),
        migrations.AddField(
            model_name='CurrentVariant',
            name='BX_ID_exLOVD',
            field=models.TextField(default=b''),
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
