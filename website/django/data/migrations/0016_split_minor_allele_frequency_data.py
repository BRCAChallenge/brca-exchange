# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations
from data.models import Variant, Report

EMPTY = '-'


def split_minor_allele_frequencies(apps, schema_editor):
    for Obj in [Variant, Report]:
        objs = Obj.objects.all()
        for obj in objs:
            if obj.Minor_allele_frequency_ESP and obj.Minor_allele_frequency_ESP != EMPTY:
                eaAlleleFrequency = EMPTY
                aaAlleleFrequency = EMPTY
                alleleFrequency = EMPTY
                maf = obj.Minor_allele_frequency_ESP.split(',')

                if len(maf) > 2:
                    alleleFrequency = "%s" % (float(maf[2]) / 100)
                if len(maf) > 1:
                    aaAlleleFrequency = "%s" % (float(maf[1]) / 100)
                if len(maf) > 0:
                    eaAlleleFrequency = "%s" % (float(maf[0]) / 100)

                obj.EA_Allele_Frequency_ESP = eaAlleleFrequency
                obj.AA_Allele_Frequency_ESP = aaAlleleFrequency
                obj.Allele_Frequency_ESP = alleleFrequency
                obj.save()


class Migration(migrations.Migration):

    dependencies = [
        ('data', '0015_auto_20170602_1458'),
    ]

    operations = [
        migrations.RunPython(split_minor_allele_frequencies),
        migrations.RunSQL(
            """
            DROP MATERIALIZED VIEW IF EXISTS currentvariant;
            CREATE MATERIALIZED VIEW currentvariant AS (
                SELECT * FROM "variant" WHERE (
                    "id" IN ( SELECT DISTINCT ON ("Genomic_Coordinate_hg38") "id" FROM "variant" ORDER BY "Genomic_Coordinate_hg38" ASC, "Data_Release_id" DESC )
                )
            );
            """
        )
    ]
