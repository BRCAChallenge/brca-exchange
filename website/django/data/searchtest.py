import json
import os
import pytest
import string
import pdb
import shutil, tempfile
from os import path
from urllib import quote
from django.http import JsonResponse, HttpResponse
from brca import settings
from data.models import Variant, CurrentVariant, ChangeType
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


TODO LIST
1. Filters (managed by filters and filter_values params)
2. Sources (managed by the include params)
3. Release (release param)
4. Change Types (change_types param)
5. Show Deleted (show_deleted param)
6. Format (format param)

Start by testing that the correct objects are used (e.g. if a release number is specified, 
the response should contain a Variant object(s), otherwise it should contain only 
CurrentVariant(s)). Next, test that each of the params mentioned in the numbered 
list above returns the expected result. To do this, you will need to make adjustments 
to existing_variant in test_data.py, and possibly create other objects with additional data.
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
        retrieved_variant = Variant.objects.get(Genomic_Coordinate_hg38=new_variant_genomic_coordinate_hg38)
        self.assertIsNotNone(retrieved_variant)
        self.assertEqual(retrieved_variant.Genomic_Coordinate_hg38, new_variant_genomic_coordinate_hg38)


    def test_current_variant_model(self):
        """This tests creation of a new Variant and that the associated CurrentVariant has the same basic properties."""
        (new_variant, new_current_variant) = create_variant_and_materialized_view(test_data.new_variant())
        new_current_variant_genomic_coordinate_hg38 = test_data.new_variant()['Genomic_Coordinate_hg38']
        retrieved_variant = CurrentVariant.objects.get(Genomic_Coordinate_hg38=new_current_variant_genomic_coordinate_hg38)
        self.assertIsNotNone(retrieved_variant)
        self.assertIsInstance(retrieved_variant, CurrentVariant)
        self.assertEqual(new_variant.Genomic_Coordinate_hg38, retrieved_variant.Genomic_Coordinate_hg38)


    def test_index_resource_json(self):
        #Tests search for all data in json format returns a JsonResponse
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

    @skip("Not complete")
    def test_format_tsv(self):
        '''Tests format parameter with format tsv'''
        #make sure there are 2 variants in database, search for 1 variant by coordinate, status code = 200, format = tsv, check response.content and test if it contains genomic coordinate and doesn't contain any other genomic coordinate
        variant_1 = test_data.existing_variant()
        variant_1['Genomic_Coordinate_hg38'] = 'non-tsv'
        (variant_1, current_variant_1) = create_variant_and_materialized_view(variant_1)
        variant_2 = test_data.existing_variant()
        variant_2['Genomic_Coordinate_hg38'] = 'tsv_format'
        (variant_2, new_variant_2) = create_variant_and_materialized_view(variant_2)

        request = self.factory.get(
            '/data/?format=tsv&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=tsv_format&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD')
        # This calls the index method from data/views.py where user queries are processed.
        response = index(request)

        self.assertIsInstance(response, HttpResponse)
        self.assertEqual(response.status_code, 200)
        #self.assertEqual(tsv.loads(response.content)['count'], 1)
        #self.assertTrue(tsv.loads(response.content)['data'][0]['Genomic_Coordinate_hg38'] != variant_2['Genomic_Coordinate_hg38'])
        self.assertEqual(response.content['count'], 1)
        self.assertTrue(response.content['data'][0]['Genomic_Coordinate_hg38'] != variant_2['Genomic_Coordinate_hg38'])

    @skip("Not complete")
    def test_format_csv(self):
        '''Tests format parameter with format csv'''
        variant_1 = test_data.existing_variant()
        variant_1['Genomic_Coordinate_hg38'] = 'non-csv'
        (variant_1, current_variant_1) = create_variant_and_materialized_view(variant_1)
        variant_2 = test_data.existing_variant()
        variant_2['Genomic_Coordinate_hg38'] = 'csv_format'
        (variant_2, new_variant_2) = create_variant_and_materialized_view(variant_2)

        request = self.factory.get(
            '/data/?format=csv&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=tsv_format&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD')
        # This calls the index method from data/views.py where user queries are processed.
        response = index(request)

        self.assertIsInstance(response, HttpResponse)
        self.assertEqual(response.status_code, 200)
        #self.assertEqual(csv.loads(response.content)['count'], 1)
        #self.assertTrue(csv.loads(response.content)['data'][0]['Genomic_Coordinate_hg38'] != variant_2['Genomic_Coordinate_hg38'])
        self.assertEqual(response.content['count'], 1)
        self.assertTrue(response.content['data'][0]['Genomic_Coordinate_hg38'] != variant_2['Genomic_Coordinate_hg38'])


    def test_search_by_id(self):
        """Tests searching for a variant by id using a filter"""
        existing_current_variant_id = self.existing_variant_materialized_view.id
        request = self.factory.get(
            '/data/?format=json&filter=id&filterValue=%s&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD' % existing_current_variant_id)
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
            '/data/?format=json&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=')
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)

        response_data = json.loads(response.content)
        self.assertEqual(response_data["count"], 0)


    def test_request_with_release_number(self):
        '''Tests that the correct objects are used when release number is specified'''
        #create a new variant with new random release ID
        new_variant_1 = test_data.new_variant()
        new_variant_1["Data_Release_id"] = 4
        (new_variant_1, new_current_variant_1) = create_variant_and_materialized_view(new_variant_1)

        #updated variant with new release ID parameter
        new_variant_2 = test_data.new_variant()
        new_variant_2["Data_Release_id"] = 5
        (new_variant_2, new_current_variant_2) = create_variant_and_materialized_view(new_variant_2)

        release_number = self.existing_variant.Data_Release_id
        request = self.factory.get(
            '/data/?format=json&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=%s&release=4&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD' % new_variant_1.Genomic_Coordinate_hg38)
        response = index(request)

        self.assertEqual(len(Variant.objects.all()), 3)
        self.assertEqual(len(CurrentVariant.objects.all()), 2)
        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)

        response_data = json.loads(response.content)
        self.assertEqual(response_data["count"], 1)

        response_variant = response_data["data"][0]

        self.assertEqual(response_variant["Data_Release_id"], new_variant_1.Data_Release_id)


    def test_request_without_release_number(self):
        '''Tests that the latest release of a variant is returned when release number is NOT specified'''     
        count = 4

        while count < 20:
            latest_new_variant = test_data.new_variant()
            latest_new_variant["Data_Release_id"] = count
            (latest_new_variant, latest_new_current_variant) = create_variant_and_materialized_view(latest_new_variant)
            count = count + 1

        request = self.factory.get(
            '/data/?format=json&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=%s&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD' % latest_new_variant.Genomic_Coordinate_hg38)
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)

        response_data = json.loads(response.content)
        self.assertEqual(response_data["count"], 1)

        response_variant = response_data["data"][0]

        self.assertEqual(response_variant["Data_Release_id"], latest_new_variant.Data_Release_id)


    #Filter testing: Gene type, pathogenicity, and source selection are all filter options
    def test_filter_by_gene_type_brca_1(self):
        '''Tests filtering by 'Gene' option by BRCA1'''
        #creates a new BRCA2 variant
        new_variant_brca2 = test_data.new_variant()
        new_variant_brca2["Gene_Symbol"] = "BRCA2"
        new_variant_brca2["Genomic_Coordinate_hg38"] = "chr13:111111:A>G"
        (new_variant_brca2, new_current_variant_brca2) = create_variant_and_materialized_view(new_variant_brca2)

        request = self.factory.get(
            '/data/?format=json&filter=Gene_Symbol&filterValue=BRCA1&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD')
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)

        response_data = json.loads(response.content)

        self.assertEqual(response_data["count"], 1)

        response_variant = response_data["data"][0]

        self.assertEqual(response_variant["Genomic_Coordinate_hg38"], self.existing_variant.Genomic_Coordinate_hg38)

    def test_filter_by_gene_type_brca_2(self):
        '''Tests filtering by 'Gene' option by BRCA2'''
        #creates a new BRCA2 variant
        new_variant_brca2 = test_data.new_variant()
        new_variant_brca2["Gene_Symbol"] = "BRCA2"
        new_variant_brca2["Genomic_Coordinate_hg38"] = "chr13:111111:A>G"
        (new_variant_brca2, new_current_variant_brca2) = create_variant_and_materialized_view(new_variant_brca2)

        request = self.factory.get(
            '/data/?format=json&filter=Gene_Symbol&filterValue=BRCA2&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD')
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)

        response_data = json.loads(response.content)
        self.assertEqual(response_data["count"], 1)

        response_variant = response_data["data"][0]
        self.assertEqual(response_variant["Genomic_Coordinate_hg38"], new_variant_brca2.Genomic_Coordinate_hg38)

    def test_filter_by_gene_type_any(self):
        '''Tests filtering by 'Gene' option by Any'''
        #creates a new BRCA2 variant
        new_variant_brca2 = test_data.new_variant()
        new_variant_brca2["Gene_Symbol"] = "BRCA2"
        new_variant_brca2["Genomic_Coordinate_hg38"] = "chr13:111111:A>G"
        (new_variant_brca2, new_current_variant_brca2) = create_variant_and_materialized_view(new_variant_brca2)

        request = self.factory.get(
            '/data/?format=json&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD')
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)

        response_data = json.loads(response.content)
        self.assertEqual(response_data["count"], 2)        


    def test_filter_by_pathogenicity(self):
        '''Tests filtering by 'Pathogenicity' option'''
        #variant with pathogenic classification
        new_variant_pathogenic = test_data.existing_variant()
        new_variant_pathogenic['Genomic_Coordinate_hg38'] = 'chr13:PATHOGENIC:A>G'
        new_variant_pathogenic['Pathogenicity_expert'] = 'Pathogenic'
        (new_variant_pathogenic, new_current_variant_pathogenic) = create_variant_and_materialized_view(new_variant_pathogenic)

        #variant with benign/little clinical significance classification
        new_variant_benign = test_data.existing_variant()
        new_variant_benign['Genomic_Coordinate_hg38'] = 'chr13:BENIGN:A>G'
        new_variant_benign['Pathogenicity_expert'] = 'Benign / Little Clinical Significance'
        (new_variant_benign, new_current_variant_benign) = create_variant_and_materialized_view(new_variant_benign)

        #variant with not_yet_reviewed status
        new_variant_not_reviewed = test_data.existing_variant()
        new_variant_not_reviewed['Genomic_Coordinate_hg38'] = 'chr13:NOT_REVIEWED:A>G'
        new_variant_not_reviewed['Pathogenicity_expert'] = 'Not Yet Reviewed'
        (new_variant_not_reviewed, new_current_variant_not_reviewed) = create_variant_and_materialized_view(new_variant_not_reviewed)

        filter_list = [
            'Pathogenic',
            'Benign / Little Clinical Significance',
            'Not Yet Reviewed'
            ]

        for filter_name in filter_list:

            message = 'Test case for Pathogenicity filter \'' + filter_name + '\' failed!!!'
            request = self.factory.get(
                '/data/?format=json&filter=Pathogenicity_expert&filterValue=%s&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD' % filter_name)
            response = index(request)

            self.assertIsInstance(response, JsonResponse, message)
            self.assertEqual(response.status_code, 200, message)

            response_data = json.loads(response.content)
            response_variants = response_data['data']
            variant_num = 1  
            for variant in response_variants:
                if self.existing_variant_materialized_view.Pathogenicity_expert == variant['Pathogenicity_expert']:
                    variant_num = 2

                self.assertEqual(variant_num, response_data['count'])
                self.assertTrue(filter_name in variant['Pathogenicity_expert'], message)


    def test_filter_by_sources(self):
        '''Tests filtering by 'Source' options'''
        new_variant_enigma = test_data.new_variant_no_source()
        new_variant_clinvar = test_data.new_variant_no_source()
        new_variant_1000_genomes = test_data.new_variant_no_source()
        new_variant_exac = test_data.new_variant_no_source()
        new_variant_lovd = test_data.new_variant_no_source()
        new_variant_bic = test_data.new_variant_no_source()
        new_variant_esp = test_data.new_variant_no_source()
        new_variant_exlovd = test_data.new_variant_no_source()

        variant_list = [
            new_variant_enigma,
            new_variant_clinvar,
            new_variant_1000_genomes,
            new_variant_exac,
            new_variant_lovd,
            new_variant_bic,
            new_variant_esp,
            new_variant_exlovd
            ]


        source_list = ['ENIGMA','ClinVar','1000_Genomes','ExAC','LOVD','BIC','ESP','exLOVD']

        '''Loops through list of new variants, assigns a new Genomic Coordinate and Source to each and 
        sets their appropriate 'Variant_in_' values to True, then puts them in database'''
        for new_variant, source in zip(variant_list, source_list):
            new_variant['Genomic_Coordinate_hg38'] = 'chr13:' + source + ':A>G'
            new_variant['Source'] = source
            new_variant['Variant_in_' + source] = True
            (new_variant, new_current_variant) = create_variant_and_materialized_view(new_variant)

        for source_name in source_list:
            message = 'Test case for source \'' + source_name + '\' failed!!!'

            request = self.factory.get(
                '/data/?format=json&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=&include=Variant_in_%s' % source_name)
            response = index(request)

            self.assertIsInstance(response, JsonResponse, message)
            self.assertEqual(response.status_code, 200, message)

            response_data = json.loads(response.content)
            #Checks if existing_variant source field is truthy
            if getattr(self.existing_variant_materialized_view, 'Variant_in_' + source_name)  == True:
                variant_num = 2
            else:
                variant_num = 1

            self.assertEqual(response_data['count'], variant_num, message)

            response_variants = response_data['data']

            for variant in response_variants:
                self.assertTrue(variant['Variant_in_' + source_name], message)


    def test_change_types(self):
        '''Tests change_types parameter'''
        change_types_list = ChangeType.objects.values()

        #creates new variants that have differing change_type_ids: one of each type
        for change_type in change_types_list:
            new_variant = test_data.new_variant()
            new_variant['Change_Type_id'] = change_type['id']
            new_variant['Data_Release_id'] = 6
            new_variant['Genomic_Coordinate_hg38'] = change_type['name']
            (new_variant, new_current_variant) = create_variant_and_materialized_view(new_variant)

        #tests that by filtering by each change_type, one (correct) result is returned
        for change_type in change_types_list:
            request = self.factory.get(
                '/data/?format=json&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=&release=6&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD&change_types=%s&show_deleted=true' % change_type['name'])
            response = index(request)

            self.assertIsInstance(response, JsonResponse)
            self.assertEqual(response.status_code, 200)

            response_data = json.loads(response.content)

            self.assertEqual(response_data['count'], 1)

            response_variant = response_data['data'][0]

            self.assertEqual(response_variant['Genomic_Coordinate_hg38'], change_type['name'])


    def test_show_deleted_not_true(self):
        '''Tests show_deleted parameter when show_deleted is not called'''
        new_variant_1 = test_data.existing_variant()
        new_variant_1['Data_Release_id'] = 7
        new_variant_1['Chr'] = 'not_deleted_1'
        (new_variant_1, new_current_variant_1) = create_variant_and_materialized_view(new_variant_1)

        new_variant_2 = test_data.existing_variant()
        new_variant_2['Data_Release_id'] = 6
        new_variant_2['Chr'] = 'not_deleted_2'
        (new_variant_2, new_current_variant_2) = create_variant_and_materialized_view(new_variant_2)

        new_variant_deleted = test_data.new_variant()
        new_variant_deleted['Data_Release_id'] = 6
        new_variant_deleted['Change_Type_id'] = 2
        new_variant_deleted['Genomic_Coordinate_hg38'] = 'deleted_variant'
        new_variant_deleted['Chr'] = 'deleted'
        (new_variant_deleted, new_current_variant_deleted) = create_variant_and_materialized_view(new_variant_deleted)

        request = self.factory.get(
            '/data/?format=json&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD')
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)

        response_data = json.loads(response.content)

        self.assertEqual(response_data['count'], 1)
        self.assertEqual(response_data['deletedCount'], 1)

        
    def test_show_deleted_true(self):
        '''Tests show_deleted parameter when show_deleted=true'''
        new_variant_1 = test_data.existing_variant()
        new_variant_1['Data_Release_id'] = 7
        new_variant_1['Chr'] = 'not_deleted_1'
        (new_variant_1, new_current_variant_1) = create_variant_and_materialized_view(new_variant_1)

        new_variant_2 = test_data.existing_variant()
        new_variant_2['Data_Release_id'] = 6
        new_variant_2['Chr'] = 'not_deleted_2'
        (new_variant_2, new_current_variant_2) = create_variant_and_materialized_view(new_variant_2)

        new_variant_deleted = test_data.new_variant()
        new_variant_deleted['Data_Release_id'] = 6
        new_variant_deleted['Change_Type_id'] = 2
        new_variant_deleted['Genomic_Coordinate_hg38'] = 'deleted_variant'
        new_variant_deleted['Chr'] = 'deleted'
        (new_variant_deleted, new_current_variant_deleted) = create_variant_and_materialized_view(new_variant_deleted)


        request = self.factory.get(
            '/data/?format=json&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD&show_deleted=true')
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)

        response_data = json.loads(response.content)

        self.assertEqual(response_data['count'], 2)
        self.assertEqual(response_data['deletedCount'], 0)


    #Begin search_by tests
    #-------------------------------------------------------------------------------------------------------------

    def test_search_by_genomic_coordinate(self):
        '''Tests that searching for a variant with a genomic_coordinate_hg38 search term is successful'''
        existing_genomic_coordinate = self.existing_variant_materialized_view.Genomic_Coordinate_hg38
        request = self.factory.get(
            '/data/?format=json&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&include=Variant_in_ENIGMA&search_term=%s' % existing_genomic_coordinate)
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)

        response_data = json.loads(response.content)
        self.assertEqual(response_data['count'], 1)

        response_variant = response_data['data'][0]
        self.assertEqual(response_variant['Genomic_Coordinate_hg38'], self.existing_variant.Genomic_Coordinate_hg38)


    def test_search_by_combined_gene_symbol_and_genomic_coordinate(self):
        '''Tests that searching for a variant with a 'gene_symbol:genomic_coordinate_hg38' search term is successful'''
        gene_symbol_genomic_coordinate_hg38 = self.existing_variant_materialized_view.Gene_Symbol + ':' + self.existing_variant_materialized_view.Genomic_Coordinate_hg38
        request = self.factory.get(
            '/data/?format=json&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=%s&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD' % gene_symbol_genomic_coordinate_hg38)
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)

        response_data = json.loads(response.content)

        self.assertEqual(response_data['count'], 1)

        response_variant = response_data['data'][0]
        self.assertEqual(response_variant['Genomic_Coordinate_hg38'], self.existing_variant.Genomic_Coordinate_hg38)


    '''
    Add tests for each of the following searches, as well as some tests to confirm that
    faulty search strings return no results.

    User submitted search --> Field:Field

    Test acceptable searches for existing variants

    Test searches for non-existing variants

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
    BRCA1:p.(Ala280Gly) --> Gene_Symbol:HGVS_Protein.split(':')[1] (HGVS_Protein is actually stored as NP_009225.1:p.(Ala280Gly), so this has to be split on the ':')
    BRCA1:A280G --> Gene_Symbol:Protein_Change
    NP_009225.1:p.(Ala280Gly) --> HGVS_Protein !!! NOT included in test_list
    NP_009225.1:A280G --> HGVS_Protein.split(':')[0]:Protein_Change
    Add test case for space separator


    EXAMPLES OF SEARCHES THAT MUST BE SUPPORTED (with space or colon delimiter)

    BRCA1:chr17:g.43094692:G>C
    BRCA1:chr17:g.41246709:G>C
    BRCA1:chr17:g.38500235:G>C
    BRCA1:958C>G
    BRCA1:c.839C>G
    NM_007294.3:chr17:g.43094692:G>C
    NM_007294.3:chr17:g.41246709:G>C
    NM_007294.3:chr17:g.38500235:G>C
    NM_007294.3:958C>G
    NM_007294.3:c.839C>G
    BRCA1:p.(Ala280Gly)
    BRCA1:A280G
    NP_009225.1:p.(Ala280Gly)
    NP_009225.1:A280G
    BRCA2 6174delT
    BRCA1 185delAG
    '''

    @skip("Full functionality not yet implemented")
    def test_search_examples(self):
        '''When search functionality is complete, this test should pass. This tests all the different ways to search'''
        test_list = [
            'BRCA1:chr17:g.43094692:G>C',
            'BRCA1:chr17:g.41246709:G>C',
            'BRCA1:chr17:g.38500235:G>C',
            'BRCA1:958C>G',
            'BRCA1:c.839C>G',
            'NM_007294.3:chr17:g.43094692:G>C',
            'NM_007294.3:chr17:g.41246709:G>C',
            'NM_007294.3:chr17:g.38500235:G>C',
            'NM_007294.3:958C>G',
            'NM_007294.3:c.839C>G',
            'BRCA1:p.(Ala280Gly)',
            'BRCA1:A280G',
            'NP_009225.1:p.(Ala280Gly)',
            'NP_009225.1:A280G',
            'BRCA2 6174delT',
            'BRCA1 185delAG'
            ]

        for test_term in test_list:
            request = self.factory.get(
            '/data/?format=json&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=%s&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD' % test_term)
            response = index(request)

            message = 'Test case for \'' + test_term + '\' failed!!!'

            self.assertIsInstance(response, JsonResponse, message)
            self.assertEqual(response.status_code, 200, message)

            response_data = json.loads(response.content)

            self.assertEqual(response_data['count'], 1, message)

            response_variant = response_data['data'][0]
            #self.assertEqual(response_variant[test_term], getattr(self.existing_variant, test_term), message)


    def test_search_by_colon_delimiters(self):
        #Tests searching for variants with colon delimiters
        test_list = {
            'Gene_Symbol':'Genomic_Coordinate_hg38',
            'Gene_Symbol':'Genomic_Coordinate_hg37',
            'Gene_Symbol':'Genomic_Coordinate_hg36',
            'Gene_Symbol':'BIC_Nomenclature',
            'Gene_Symbol':'HGVS_cDNA',
            'Reference_Sequence':'Genomic_Coordinate_hg38',
            'Reference_Sequence':'Genomic_Coordinate_hg37',
            'Reference_Sequence':'Genomic_Coordinate_hg36',
            'Reference_Sequence':'BIC_Nomenclature',
            'Reference_Sequence':'HGVS_cDNA',
            'Gene_Symbol':'Protein_Change'
        }

        #This loops through test_list dictionary, providing matching pairs of correctly-ordered search terms as field_1 and field_2
        for field_1, field_2 in test_list.items():
            test_case = getattr(self.existing_variant_materialized_view, field_1) + ':' + getattr(self.existing_variant_materialized_view, field_2)

            request = self.factory.get(
            '/data/?format=json&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=%s&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD' % test_case)
            response = index(request)
            #build a string for failure message
            message = 'Test case for \'' + field_1 + ':' + field_2 + '\' failed!!!'

            self.assertIsInstance(response, JsonResponse, message)
            self.assertEqual(response.status_code, 200, message)

            response_data = json.loads(response.content)

            self.assertEqual(response_data['count'], 1, message)

            response_variant = response_data['data'][0]
            self.assertEqual(response_variant[field_2], getattr(self.existing_variant,field_2), message)


    def test_search_by_gene_symbol_and_hgvs_protein_with_colon(self):
        '''Tests searching for 'BRCA1 p.(Ala280Gly) --> Gene_Symbol HGVS_Protein.split(':')[1]' with colon separator
        BRCA1:p.(Ala280Gly) --> Gene_Symbol:HGVS_Protein.split(':')[1] (HGVS_Protein is actually stored as NP_009225.1:p.(Ala280Gly), so this has to be split on the ':')'''
        gene_symbol_hgvs_protein = self.existing_variant_materialized_view.Gene_Symbol + ':' + self.existing_variant_materialized_view.HGVS_Protein.split(':')[1]
        request = self.factory.get(
            '/data/?format=json&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=%s&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD' % gene_symbol_hgvs_protein)
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)

        response_data = json.loads(response.content)

        self.assertEqual(response_data['count'], 1)

        response_variant = response_data['data'][0]
        self.assertEqual(response_variant['HGVS_Protein'], self.existing_variant.HGVS_Protein)


    def test_search_by_hgvs_protein_and_protein_change_with_colon(self):
        '''Tests searching for 'HGVS_Protein:Protein_Change' with colon separator
        NP_009225.1:A280G --> HGVS_Protein.split(':')[0] Protein_Change'''
        hgvs_protein_protein_change = self.existing_variant_materialized_view.HGVS_Protein.split(':')[0] + ':' + self.existing_variant_materialized_view.Protein_Change
        request = self.factory.get(
            '/data/?format=json&order_by=HGVS_Protein&direction=ascending&page_size=20&page_num=0&search_term=%s&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD' % hgvs_protein_protein_change)
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)

        response_data = json.loads(response.content)

        self.assertEqual(response_data['count'], 1)

        response_variant = response_data['data'][0]
        self.assertEqual(response_variant['Protein_Change'], self.existing_variant.Protein_Change)


    @skip('Spaces not yet implemented')
    def test_search_by_space_delimiters(self):
        #Tests searching for different search terms with space delimiters
        test_list = {
            'Gene_Symbol':'Genomic_Coordinate_hg38',
            'Gene_Symbol':'Genomic_Coordinate_hg37',
            'Gene_Symbol':'Genomic_Coordinate_hg36',
            'Gene_Symbol':'BIC_Nomenclature',
            'Gene_Symbol':'HGVS_cDNA',
            'Reference_Sequence':'Genomic_Coordinate_hg38',
            'Reference_Sequence':'Genomic_Coordinate_hg37',
            'Reference_Sequence':'Genomic_Coordinate_hg36',
            'Reference_Sequence':'BIC_Nomenclature',
            'Reference_Sequence':'HGVS_cDNA',
            'Gene_Symbol':'Protein_Change'
        }
        #This loops through the test_list dictionary, providing matching pairs of correctly-ordered search terms
        for field_1, field_2 in test_list.items():
            test_case = getattr(self.existing_variant_materialized_view, field_1) + ' ' + getattr(self.existing_variant_materialized_view, field_2)

            request = self.factory.get(
            '/data/?format=json&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=%s&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD' % test_case)
            response = index(request)
            #build a string for failure message
            message = 'Test case for \'' + field_1 + ' ' + field_2 + '\' failed!!!'

            self.assertIsInstance(response, JsonResponse, message)
            self.assertEqual(response.status_code, 200, message)

            response_data = json.loads(response.content)

            self.assertEqual(response_data['count'], 1, message)

            response_variant = response_data['data'][0]
            self.assertEqual(response_variant[field_2], getattr(self.existing_variant,field_2), message)


    @skip('Spaces not yet implemented')
    def test_search_by_gene_symbol_and_hgvs_protein_with_space(self):
        '''Tests searching for 'BRCA1 p.(Ala280Gly) --> Gene_Symbol HGVS_Protein.split(':')[1]' with colon separator
        BRCA1:p.(Ala280Gly) --> Gene_Symbol:HGVS_Protein.split(':')[1] (HGVS_Protein is actually stored as NP_009225.1:p.(Ala280Gly), so this has to be split on the ':')'''
        gene_symbol_hgvs_protein = self.existing_variant_materialized_view.Gene_Symbol + ' ' + self.existing_variant_materialized_view.HGVS_Protein.split(':')[1]
        request = self.factory.get(
            '/data/?format=json&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=%s&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD' % gene_symbol_hgvs_protein)
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)

        response_data = json.loads(response.content)

        self.assertEqual(response_data['count'], 1)

        response_variant = response_data['data'][0]
        self.assertEqual(response_variant['HGVS_Protein'], self.existing_variant.HGVS_Protein)


    @skip('Spaces not yet implemented')
    def test_search_by_hgvs_protein_and_protein_change_with_space(self):
        '''Tests searching for 'HGVS_Protein:Protein_Change' with colon separator
        NP_009225.1:A280G --> HGVS_Protein.split(':')[0] Protein_Change'''
        hgvs_protein_protein_change = self.existing_variant_materialized_view.HGVS_Protein.split(':')[0] + ':' + self.existing_variant_materialized_view.Protein_Change
        request = self.factory.get(
            '/data/?format=json&order_by=HGVS_Protein&direction=ascending&page_size=20&page_num=0&search_term=%s&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD' % hgvs_protein_protein_change)
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)

        response_data = json.loads(response.content)

        self.assertEqual(response_data['count'], 1)

        response_variant = response_data['data'][0]
        self.assertEqual(response_variant['Protein_Change'], self.existing_variant.Protein_Change)


    @skip('Functionality not yet implemented')
    def test_search_by_combined_gene_symbol_and_genomic_coordinate_with_space(self):
        '''Tests searching for 'gene_symbol genomic_coordinate_hg38' is successful

        NOTE: This is an example of a test for future functionality. It fails (as it should) because we have not implemented code
        to support this search type (only colon delimiters at the moment, not spaces). Failing tests are incredibly useful
        to developers when they implement a new feature that would cause the test to pass. Writing tests to define the
        parameters of a new feature prior to implementation is known as Test Driven Development (TDD).

        If you find the failure bothersome, uncomment the skip line above the test.'''
        gene_symbol_genomic_coordinate_hg38 = self.existing_variant_materialized_view.Gene_Symbol + ' ' + self.existing_variant_materialized_view.Genomic_Coordinate_hg38
        request = self.factory.get(
            '/data/?format=json&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=%s&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD' % gene_symbol_genomic_coordinate_hg38)
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)

        response_data = json.loads(response.content)
        self.assertEqual(response_data['count'], 1)

        response_variant = response_data['data'][0]
        self.assertEqual(response_variant['Genomic_Coordinate_hg38'], self.existing_variant.Genomic_Coordinate_hg38)



    #FAULTY STRING SEARCHES (These should fail/return no result)
    #**********************************************************************************
    def test_search_by_incorrect_delimiters(self):
        '''Tests searches that incorrect delimiters do not work'''
        symbol_list = [';','.','+','-','>','(',')','\'']

        for delimiter in symbol_list:
            gene_symbol_genomic_coordinate_hg38 = self.existing_variant_materialized_view.Gene_Symbol + delimiter + self.existing_variant_materialized_view.Genomic_Coordinate_hg38
            request = self.factory.get(
                '/data/?format=json&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=%s' % gene_symbol_genomic_coordinate_hg38)
            response = index(request)

            message = 'Delimiter test case with \'' + delimiter + '\' fails!!!'

            self.assertIsInstance(response, JsonResponse, message)
            self.assertEqual(response.status_code, 200, message)

            response_data = json.loads(response.content)

            self.assertEqual(response_data['count'], 0, message)


    def test_search_by_gibberish(self):
        '''Tests gibberish searches
        gibberish,<.>/?'';:[{]}\|=+-_)(*&%^$#@!~`'''
        ascii_list = [chr(i) for i in xrange(256)]
        ascii_string = 'hey I\'m gibberish look at meeee' + ''.join(ascii_list)

        request = self.factory.get(
            '/data/?format=json&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=%s' % ascii_string)
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)

        response_data = json.loads(response.content)
        self.assertEqual(response_data['count'], 0)
        

    def test_search_by_too_many_fields(self):
        '''Tests searches that input more than two fields
        Field:Field:Field'''
        hgvs_protein_protein_change_gene_symbol = ':' + self.existing_variant_materialized_view.HGVS_Protein.split(':')[0] + ':' + self.existing_variant_materialized_view.Protein_Change + ':' + self.existing_variant_materialized_view.Gene_Symbol
        request = self.factory.get(
            '/data/?format=json&order_by=HGVS_Protein&direction=ascending&page_size=20&page_num=0&search_term=%s&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD' % hgvs_protein_protein_change_gene_symbol)
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)

        response_data = json.loads(response.content)

        self.assertEqual(response_data['count'], 0)


    def test_search_by_too_many_colons(self):
        '''Tests searches with more than one colon in a row
        Field::Field'''
        hgvs_protein_protein_change = self.existing_variant_materialized_view.HGVS_Protein.split(':')[0] + ':::' + self.existing_variant_materialized_view.Protein_Change
        request = self.factory.get(
            '/data/?format=json&order_by=HGVS_Protein&direction=ascending&page_size=20&page_num=0&search_term=%s&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD' % hgvs_protein_protein_change)
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)

        response_data = json.loads(response.content)

        self.assertEqual(response_data['count'], 0)


    def test_search_by_beginning_colon(self):
        '''Tests searches with a beginning colon (uses hgvs_protein and protein_change fields as a test case)
        :Field:Field'''
        hgvs_protein_protein_change = ':' + self.existing_variant_materialized_view.HGVS_Protein.split(':')[0] + ':' + self.existing_variant_materialized_view.Protein_Change
        request = self.factory.get(
            '/data/?format=json&order_by=HGVS_Protein&direction=ascending&page_size=20&page_num=0&search_term=%s&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD' % hgvs_protein_protein_change)
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)

        response_data = json.loads(response.content)

        self.assertEqual(response_data['count'], 0)


    @skip('Not complete')
    def test_search_by_sql_injection(self):
        '''Tests the safeguards against SQL Injection
        '''

    def test_search_by_flipped_search_terms(self):
        #tests that flipping search terms does not work
        test_list = {
            'Gene_Symbol':'Genomic_Coordinate_hg38',
            'Gene_Symbol':'Genomic_Coordinate_hg37',
            'Gene_Symbol':'Genomic_Coordinate_hg36',
            'Gene_Symbol':'BIC_Nomenclature',
            'Gene_Symbol':'HGVS_cDNA',
            'Reference_Sequence':'Genomic_Coordinate_hg38',
            'Reference_Sequence':'Genomic_Coordinate_hg37',
            'Reference_Sequence':'Genomic_Coordinate_hg36',
            'Reference_Sequence':'BIC_Nomenclature',
            'Reference_Sequence':'HGVS_cDNA',  
            'Gene_Symbol':'Protein_Change'
        }

        for field_1, field_2 in test_list.items():
            test_case = getattr(self.existing_variant_materialized_view,field_2) + ':' + getattr(self.existing_variant_materialized_view,field_1)

            request = self.factory.get(
            '/data/?format=json&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=%s&include=Variant_in_ENIGMA&include=Variant_in_ClinVar&include=Variant_in_1000_Genomes&include=Variant_in_ExAC&include=Variant_in_LOVD&include=Variant_in_BIC&include=Variant_in_ESP&include=Variant_in_exLOVD' % test_case)
            response = index(request)

            #build a string for failure message
            message = 'Test case for (flipped) ' + field_2 + ':' + field_1 + ' failed!!!'

            self.assertIsInstance(response, JsonResponse, message)
            self.assertEqual(response.status_code, 200, message)

            response_data = json.loads(response.content)

            self.assertEqual(response_data['count'], 0, message)


