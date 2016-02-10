from __future__ import unicode_literals

from django.db import models


class VariantManager(models.Manager):
    use_in_migrations = True

    def create_variant(self, row):
        return self.create(gene_symbol=row[0],
                           genomic_coordinate=row[1],
                           reference_sequence=row[2],
                           hgvs_cdna=row[3],
                           bic_nomenclature=row[4],
                           abbrev_aa_change=row[5],
                           url=row[6],
                           condition_id_type=row[7],
                           condition_id_value=row[8],
                           condition_category=row[9],
                           clinical_significance=row[10],
                           date_last_evaluated=row[11],
                           assertion_method=row[12],
                           assertion_method_citation=row[13],
                           clinical_significance_citations=row[14],
                           comment_on_clinical_significance=row[15],
                           collection_method=row[16],
                           allele_origin=row[17],
                           clinvaraccession=row[18],
                           hgvs_protein=row[19])


class Variant(models.Model):
    gene_symbol = models.CharField(max_length=512)
    genomic_coordinate = models.CharField(max_length=512)
    reference_sequence = models.CharField(max_length=512)
    hgvs_cdna = models.CharField(max_length=512)
    bic_nomenclature = models.CharField(max_length=512)
    abbrev_aa_change = models.CharField(max_length=512)
    url = models.CharField(max_length=512)
    condition_id_type = models.CharField(max_length=512)
    condition_id_value = models.CharField(max_length=512)
    condition_category = models.CharField(max_length=512)
    clinical_significance = models.CharField(max_length=512)
    date_last_evaluated = models.CharField(max_length=512)
    assertion_method = models.CharField(max_length=512)
    assertion_method_citation = models.CharField(max_length=512)
    clinical_significance_citations = models.CharField(max_length=512)
    comment_on_clinical_significance = models.CharField(max_length=512)
    collection_method = models.CharField(max_length=512)
    allele_origin = models.CharField(max_length=512)
    clinvaraccession = models.CharField(max_length=512)
    hgvs_protein = models.CharField(max_length=512)

    objects = VariantManager()
