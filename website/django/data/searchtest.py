import json
import os
import pytest
from urllib import quote
from django.http import JsonResponse, HttpResponse
from brca import settings
from data.models import Variant, CurrentVariant
from django.test import TestCase, RequestFactory
import data.views as views
from django.db import connection
from unittest import skip
from django.test.client import RequestFactory
from data import test_data
from data.views import index, autocomplete


'''
NOTE:
We use materialized views of variants to handle user queries more efficiently.
Because of this, we need to test both the Variant and CurrentVariant models.
The CurrentVariant model is a materialized view of the Variant model.
see https://www.postgresql.org/docs/9.3/static/rules-materializedviews.html for
more information about materialized views. The models can be found in data/models.py.
'''


'''
TODO:
Fix the broken tests in this file. To fix the tests, DO NOT CHANGE the assertions, instead
change the code prior to the assertions so the assertions are truthy.
'''


def create_variant_and_materialized_view(variant_data):
    """
    This is a convenience method that creates/saves a new Variant and refreshes the materialized view
    at the same time (meaning a new CurrentVariant is created to match the Variant).
    """
    variant = Variant.objects.create_variant(row=variant_data)
    with connection.cursor() as cursor:
        cursor.execute("REFRESH MATERIALIZED VIEW currentvariant")
    materialized_view = CurrentVariant.objects.get(Genomic_Coordinate_hg38=variant.Genomic_Coordinate_hg38)
    return (variant, materialized_view)


class VariantTestCase(TestCase):
    def setUp(self):
        """
        All tests have access to anything created in this method. Right now, all we need is to set up a request factory
        and create a variant and a materialized view of that same variant. The variant data is pulled from the
        test_data.py file in this same directory.
        """
        self.factory = RequestFactory()
        (self.existing_variant, self.existing_variant_materialized_view) = create_variant_and_materialized_view(test_data.existing_variant())

    def test_variant_model(self):
        """This tests creation and retreival of a new variant by the Genomic_Coordinate_hg38 column"""
        Variant.objects.create_variant(row=(test_data.new_variant()))
        new_variant_genomic_coordinate_hg38 = test_data.new_variant()['Genomic_Coordinate_hg38']
        retrieved_variant = Variant.objects.get(Genomic_Coordinate_hg38='chr17:g.43070959:A>G')
        self.assertIsNotNone(retrieved_variant)
        self.assertEqual(retrieved_variant.Genomic_Coordinate_hg38, new_variant_genomic_coordinate_hg38)

    def test_current_variant_model(self):
        """This tests creation of a new Variant and that the associated CurrentVariant has the same basic properties."""
        (new_variant, new_current_variant) = create_variant_and_materialized_view(test_data.new_variant())
        new_current_variant_genomic_coordinate_hg38 = test_data.new_variant()['Genomic_Coordinate_hg38']
        retrieved_variant = Variant.objects.get(Genomic_Coordinate_hg38=new_current_variant_genomic_coordinate_hg38)
        self.assertIsNotNone(retrieved_variant)
        self.assertIsInstance(retrieved_variant, CurrentVariant)
        self.assertEqual(new_variant.Genomic_Coordinate_hg38, retrieved_variant.Genomic_Coordinate_hg38)

    def test_index_resource_json(self):
        """
        Tests search for all data in json format returns a JsonResponse
        Hint: how does the index method in views.py handle Variant models vs. CurrentVariant models?
        Double hint: what type of model is self.existing_variant?
        """
        request = self.factory.get(
            '/data/?format=json&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD')
        request.variant = self.existing_variant_materialized_view

        # This calls the index method from data/views.py where user queries are processed.
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)
        self.assertEqual(json.loads(response.content)['count'], 1)

    def test_index_resource_csv(self):
        """Tests searching for all the data in csv format returns an HttpResponse"""
        request = self.factory.get(
            '/data/?format=json&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD')
        response = index(request)

        self.assertIsInstance(response, HttpResponse)
        self.assertEqual(response.status_code, 200)

    def test_search_by_id(self):
        """Tests searching for a variant by id using a filter"""
        existing_current_variant_id = self.existing_variant_materialized_view.id
        request = self.factory.get(
            '/data/?format=json&filter=id&filterValue=%s&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD' % 5)
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)

        response_data = json.loads(response.content)
        self.assertEqual(response_data["count"], 1)

        response_variant = response_data["data"][0]
        self.assertEqual(response_variant["Genomic_Coordinate_hg38"], self.existing_variant.Genomic_Coordinate_hg38)

    def test_source_filters_all_off(self):
        """Tests all source filters on returns no variants"""
        request = self.factory.get(
            '/data/?format=json&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=')
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)

        response_data = json.loads(response.content)
        self.assertEqual(response_data["count"], 0)

    def test_source_filter_for_variant_excludes_variant(self):
        """Tests all source filters on returns no variants"""
        request = self.factory.get(
            '/data/?format=json&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=&include=Variant_in_ENIGMA')
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)

        response_data = json.loads(response.content)
        self.assertEqual(response_data["count"], 0)

    def test_search_by_genomic_coordinate(self):
        '''Tests that searching for a variant with a genomic_coordinate_hg38 search term is successful'''
        existing_genomic_coordinate = self.existing_variant_materialized_view.Genomic_Coordinate_hg38
        request = self.factory.get(
            '/data/?format=json&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&include=Variant_in_ENIGMA&search_term=%s' % "chr17:999999:A>G")
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)

        response_data = json.loads(response.content)
        self.assertEqual(response_data["count"], 1)

        response_variant = response_data["data"][0]
        self.assertEqual(response_variant["Genomic_Coordinate_hg38"], self.existing_variant.Genomic_Coordinate_hg38)

    def test_search_by_combined_gene_symbol_and_genomic_coordinate(self):
        """Tests that searching for a variant with a 'gene_symbol:genomic_coordinate_hg38' search term is successful"""
        gene_symbol_genomic_coordinate_hg38 = self.existing_variant_materialized_view.Gene_Symbol + '_' + self.existing_variant_materialized_view.Genomic_Coordinate_hg38
        request = self.factory.get(
            '/data/?format=json&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=%s&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD' % gene_symbol_genomic_coordinate_hg38)
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)

        response_data = json.loads(response.content)
        self.assertEqual(response_data["count"], 1)

        response_variant = response_data["data"][0]
        self.assertEqual(response_variant["Genomic_Coordinate_hg38"], self.existing_variant.Genomic_Coordinate_hg38)

    @skip("Functionality not yet implemented")
    def test_search_by_combined_gene_symbol_and_genomic_coordinate_with_space(self):
        """Tests searching for 'gene_symbol genomic_coordinate_hg38' is successful"""
        """
        NOTE: This is an example of a test for future functionality. It fails (as it should) because we have not implemented code
        to support this search type (only colon delimiters at the moment, not spaces). Failing tests are incredibly useful
        to developers when they implement a new feature that would cause the test to pass. Writing tests to define the
        parameters of a new feature prior to implementation is known as Test Driven Development (TDD).

        If you find the failure bothersome, uncomment the skip line above the test.
        """
        gene_symbol_genomic_coordinate_hg38 = self.existing_variant_materialized_view.Gene_Symbol + ' ' + self.existing_variant_materialized_view.Genomic_Coordinate_hg38
        request = self.factory.get(
            '/data/?format=json&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=%s&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD' % gene_symbol_genomic_coordinate_hg38)
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)

        response_data = json.loads(response.content)
        self.assertEqual(response_data["count"], 1)

        response_variant = response_data["data"][0]
        self.assertEqual(response_variant["Genomic_Coordinate_hg38"], self.existing_variant.Genomic_Coordinate_hg38)

    '''
    TODO:

    Add tests for each of the following searches, as well as some tests to confirm that
    faulty search strings return no results.

    User submitted search --> Field:Field

    BRCA1:chr17:g.43094692:G>C --> Gene_Symbol:Genomic_Coordinate_hg38
    BRCA1:chr17:g.41246709:G>C --> Gene_Symbol:Genomic_Coordinate_hg37
    BRCA1:chr17:g.38500235:G>C --> Gene_Symbol:Genomic_Coordinate_hg36
    BRCA1:958C>G --> Gene_Symbol:BIC_Nomenclature
    BRCA1:c.839C>G --> Gene_Symbol:HGVS_cDNA
    NM_007294.3:chr17:g.43094692:G>C --> Reference_Sequence:Genomic_Coordinate_hg38
    NM_007294.3:chr17:g.41246709:G>C --> Reference_Sequence:Genomic_Coordinate_hg37
    NM_007294.3:chr17:g.38500235:G>C --> Reference_Sequence:Genomic_Coordinate_hg36
    NM_007294.3:958C>G --> Reference_Sequence:BIC_Nomenclature
    NM_007294.3:c.839C>G --> Reference_Sequence:HGVS_cDNA
    BRCA1:p.(Ala280Gly) --> Gene_Symbol:HGVS_Protein.split(':')[1] (HGVS_Protein is actually stored as NP_009225.1:p.(Ala280Gly), so this has to be split on the ":")
    BRCA1:A280G --> Gene_Symbol:Protein_Change
    NP_009225.1:p.(Ala280Gly) --> HGVS_Protein
    NP_009225.1:A280G --> HGVS_Protein.split(':')[0]:Protein_Change
    '''
