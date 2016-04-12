from __future__ import unicode_literals

from django.db import migrations


class Migration(migrations.Migration):
    dependencies = [
        ('data', '0003_data'),
    ]

    operations = [
        migrations.RunSQL("""
        CREATE TABLE words AS SELECT DISTINCT left(word, 300) as word FROM (
        SELECT regexp_split_to_table(lower("Genomic_Coordinate_hg38"), '[\s|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("Genomic_Coordinate_hg37"), '[\s|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("Genomic_Coordinate_hg36"), '[\s|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("Condition_category_ENIGMA"), '[\s|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("Clinical_significance_ENIGMA"), '[\s|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("Gene_Symbol"), '[\s|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("Reference_Sequence"), '[\s|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("HGVS_cDNA"), '[\s|:''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("BIC_Identifier"), '[\s|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("HGVS_Protein"), '[\s|''"]') as word from variant
        )
        AS combined_words;

        CREATE INDEX words_idx ON words(word text_pattern_ops);
    """)

    ]
