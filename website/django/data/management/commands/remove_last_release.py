from django.core.management.base import BaseCommand, CommandError
from django.db import connection
from data.models import Variant, DataRelease


class Command(BaseCommand):
    help = 'Remove the most recent release from the database'

    def update_autocomplete_words(self):
        # Drop words table and recreate with latest data
        with connection.cursor() as cursor:
            cursor.execute("""
                DROP TABLE IF EXISTS words;
                CREATE TABLE words AS SELECT DISTINCT left(word, 300) as word, release_id FROM (
                SELECT regexp_split_to_table(lower("Genomic_Coordinate_hg38"), '[\s|''"]') as word, "Data_Release_id" as release_id from variant UNION
                SELECT regexp_split_to_table(lower("Genomic_Coordinate_hg37"), '[\s|''"]') as word, "Data_Release_id" as release_id from variant UNION
                SELECT regexp_split_to_table(lower("Genomic_Coordinate_hg36"), '[\s|''"]') as word, "Data_Release_id" as release_id from variant UNION
                SELECT regexp_split_to_table(lower("Clinical_significance_ENIGMA"), '[\s|''"]') as word, "Data_Release_id" as release_id from variant UNION
                SELECT regexp_split_to_table(lower("Gene_Symbol"), '[\s|''"]') as word, "Data_Release_id" as release_id from variant UNION
                SELECT regexp_split_to_table(lower("Reference_Sequence"), '[\s|''"]') as word, "Data_Release_id" as release_id from variant UNION
                SELECT regexp_split_to_table(lower("HGVS_cDNA"), '[\s|:''"]')  as word, "Data_Release_id" as release_id from variant UNION
                SELECT regexp_split_to_table(lower("BIC_Nomenclature"), '[\s|''"]') as word, "Data_Release_id" as release_id from variant UNION
                SELECT regexp_split_to_table(lower("HGVS_Protein"), '[\s|''"]') as word, "Data_Release_id" as release_id from variant
                )
                AS combined_words;

                CREATE INDEX words_idx ON words(word text_pattern_ops);
            """)

    def handle(self, *args, **options):
        latest_release_id = DataRelease.objects.all().order_by("-id")[0]

        Variant.objects.filter(Data_Release_id=latest_release_id).delete()
        DataRelease.objects.filter(id=latest_release_id).delete()

        self.update_autocomplete_words()

        # update materialized view of current variants
        with connection.cursor() as cursor:
            cursor.execute("REFRESH MATERIALIZED VIEW currentvariant")
