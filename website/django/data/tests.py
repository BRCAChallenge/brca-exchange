import json
import os
import unittest
from urllib import quote

from django.http import JsonResponse, HttpResponse
from django.test import TestCase
from django.test.client import RequestFactory

from brca import settings
from data import test_data
from data.models import Variant
from data.views import index, autocomplete


class VariantTestCase(TestCase):
    def setUp(self):
        self.factory = RequestFactory()
        datafile = os.path.join(settings.BASE_DIR, 'data', 'resources', 'aggregated.tsv')
        self.db_size = sum(1 for _ in open(datafile)) - 1

    def test_variant_model(self):
        """Create a new variant and then retrieve it by the Genomic_Coordinate_hg38 column"""
        self.assertEqual(len(Variant.objects.all()), self.db_size)
        Variant.objects.create_variant(row=(test_data.new_variant()))
        self.assertEqual(len(Variant.objects.all()), self.db_size + 1)
        retrieved_variant = Variant.objects.get(Genomic_Coordinate_hg38="chr17:999999:A>G")
        self.assertIsNotNone(retrieved_variant)

    def test_index_resource_json(self):
        """Searching for all the data in json format returns a JsonResponse"""
        request = self.factory.get(
            '/data/?format=json&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=')
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)
        self.assertEqual(json.loads(response.content)['count'], self.db_size)

    def test_index_resource_csv(self):
        """Searching for all the data in csv format returns an HttpResponse"""
        request = self.factory.get(
            '/data/?format=json&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=')
        response = index(request)

        self.assertIsInstance(response, HttpResponse)
        self.assertEqual(response.status_code, 200)

    def search_by_id(self):
        """Searching for a variant by id using a filter should return the expected result"""
        existing_id = test_data.existing_variant()["id"]
        request = self.factory.get(
            '/data/?format=json&filter=id&filterValue=%s&order_by=Gene_Symbol&direction=ascending&page_size=20&page_num=0&search_term=' % existing_id)
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)
        self.assertJSONEqual(response.content, {"count": 1, "data": [test_data.existing_variant()]})

    def test_autocomplete_nuclotide(self):
        """Getting autocomplete suggestions for words starting with c.2123 should return 2 results"""
        search_term = quote('c.2123')
        expected_autocomplete_results = [["c.2123c>a"], ["c.2123c>t"]]

        request = self.factory.get('/data/suggestions/?term=%s' % search_term)
        response = autocomplete(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)
        self.assertJSONEqual(response.content, {"suggestions": expected_autocomplete_results})

    def test_autocomplete_bic(self):
        """Getting autocomplete suggestions for words starting with IVS7+10 should return 2 results"""
        search_term = quote('ivs7+10')
        expected_autocomplete_results = [["ivs7+1028t>a"], ["ivs7+1037t>c"]]

        query = '/data/suggestions/?term=%s' % search_term
        request = self.factory.get(query)
        response = autocomplete(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)
        self.assertJSONEqual(response.content, {"suggestions": expected_autocomplete_results})


if __name__ == '__main__':
    unittest.main()
