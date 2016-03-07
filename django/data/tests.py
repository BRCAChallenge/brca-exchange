import json
import unittest

from django.http import JsonResponse, HttpResponse
from django.test import TestCase
from django.test.client import RequestFactory

from data import test_data
from data.models import Variant
from data.views import index, autocomplete


class VariantTestCase(TestCase):
    def setUp(self):
        self.factory = RequestFactory()

    def test_variant_model(self):
        """Create a new variant and then retrieve it by HGVS_genomic"""
        self.assertEqual(len(Variant.objects.all()), 13447)
        Variant.objects.create_variant(row=(test_data.new_variant()))
        self.assertEqual(len(Variant.objects.all()), 13448)
        retrieved_variant = Variant.objects.get(HGVS_genomic="NC_000017.10:g.1234567T>A")
        self.assertIsNotNone(retrieved_variant)

    def test_index_resource_json(self):
        """Searching for all the data in json format returns a JsonResponse"""
        request = self.factory.get(
            '/data/?format=json&order_by=Gene_symbol&direction=ascending&page_size=20&page_num=0&search_term=')
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)
        self.assertEqual(json.loads(response.content)['count'], 13447)

    def test_index_resource_csv(self):
        """Searching for all the data in csv format returns an HttpResponse"""
        request = self.factory.get(
            '/data/?format=csv&order_by=Gene_symbol&direction=ascending&page_size=20&page_num=0&search_term=')
        response = index(request)

        self.assertIsInstance(response, HttpResponse)
        self.assertEqual(response.status_code, 200)

    def search_by_id(self):
        """Searching for a variant by id using a filter should return the expected result"""
        request = self.factory.get(
            '/data/?format=json&filter=id&filterValue=5456&order_by=Gene_symbol&direction=ascending&page_size=1&page_num=0&search_term=')
        response = index(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)
        self.assertJSONEqual(response.content, {"count": 1, "data": [test_data.existing_variant()]})

    def test_autocomplete(self):
        """Getting autocomplete suggestions for words starting with c.2123 should return 2 results"""
        request = self.factory.get('/data/suggestions/?term=c.2123')
        response = autocomplete(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)
        self.assertJSONEqual(response.content, {"suggestions": [["c.2123c>a"], ["c.2123c>t"]]})


if __name__ == '__main__':
    unittest.main()
