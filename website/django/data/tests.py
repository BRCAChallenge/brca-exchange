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

##################### MY EDITS #########################
from data.views import variant_search, brca_to_ga4gh, ErrorMessages, get_offset, get_var_by_id, get_variantSet, get_varset_by_id, varsetId_empty_catcher, empty_varId_catcher, search_datasets
from django.test import Client
c = Client()
import google.protobuf.json_format as json_format
from ga4gh import variant_service_pb2 as v_s
from ga4gh import variants_pb2 as vrs

##################### MY EDITS END#####################

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

    #@unittest.skip("Not Passing")
    def test_autocomplete_nucleotide(self):
        """Getting autocomplete suggestions for words starting with c.2123 should return 2 results"""
        search_term = quote('c.2123')
        expected_autocomplete_results = [["c.2123c>a"], ["c.2123c>t"]]

        request = self.factory.get('/data/suggestions/?term=%s' % search_term)
        response = autocomplete(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)
        self.assertJSONEqual(response.content, {"suggestions": expected_autocomplete_results})
    #@unittest.skip("Not Passing")
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

###################################################################################
############################# NEW TESTS START #####################################
    def test_ga4gh_variants_status_code(self):
        request0 =  self.factory.post(
            "/data/ga4gh/variants/search", json.dumps({"referenceName": "chr17", "variantSetId": "brca-hg37", "start": 51425158, "end": 51425158029, "pageSize": 5 }), content_type="application/json")
        response = variant_search(request0)
        self.assertEqual(response.status_code, 200)

    def test_ga4gh_translator(self):
        var_resp = vrs.Variant()
        resp = Variant.objects.values()[0]
        genRefer = "hg37"
        for j in resp:
            if resp[j] == "-" or (resp[j] == ""):
                continue
            if j == "Genomic_Coordinate_"+genRefer:
                var_resp.reference_name, start, bases = resp[j].split(':')
                var_resp.reference_bases, alternbases = bases.split(">")
                for i in range(len(alternbases)):
                    var_resp.alternate_bases.append(alternbases[i])
                var_resp.start = int(start)
                var_resp.end = var_resp.start + len(alternbases)
                continue
            if j == "id":
                var_resp.id = genRefer+"-"+str(resp['id'])
                var_resp.variant_set_id = "brca-exchange"+"-"+genRefer

                var_resp.created = 0
                var_resp.updated = 0
                continue
            if j == "Synonyms":
                Names = [i for i in str(resp[j]).split(",")]
                for i in Names:
                    var_resp.names.append(i)
            else:
                var_resp.info[str(j)].append(resp[j])

        expectedResp = json.dumps(json_format._MessageToJsonObject(var_resp, False))
        jresponse = json.dumps(json_format._MessageToJsonObject(brca_to_ga4gh(resp, genRefer), False))
        self.assertJSONEqual(jresponse , expectedResp)

    def test_validated_request(self):
        request = v_s.SearchVariantsRequest()
        req = self.factory.post("/data/ga4gh/variants/search", json.dumps(json_format._MessageToJsonObject(request, False)),
                                content_type="application/json")
        response = variant_search(req)
        self.assertJSONEqual(response.content, ErrorMessages['VariantSetId'])

        request.variant_set_id = "Something not null"
        req = self.factory.post("/data/ga4gh/variants/search", json.dumps(json_format._MessageToJsonObject(request, False)),
                                content_type="application/json")
        response = variant_search(req)
        self.assertJSONEqual(response.content, ErrorMessages['referenceName'])

        request.reference_name = "chr17"
        req= self.factory.post("/data/ga4gh/variants/search", json.dumps(json_format._MessageToJsonObject(request, False)),
                                content_type="application/json")
        response = variant_search(req)
        self.assertJSONEqual(response.content, ErrorMessages['start'])

        request.start = 14589
        req = self.factory.post("/data/ga4gh/variants/search", json.dumps(json_format._MessageToJsonObject(request, False)),
                                content_type="application/json")
        response = variant_search(req)
        self.assertJSONEqual(response.content, ErrorMessages['end'])

    def test_offset_calculator(self):
        start = 12500
        end = 13000
        s1, e1 = get_offset(start, end) #Ending offset is the same, the query will be exclusive, and therefore query will satisfy type, that is inclusive.
        self.assertEquals(s1, start+1)
        self.assertEquals(e1, end+1)
        opperlen = ['a','b']
        s2, e2 = get_offset(start, end, opperlen)
        self.assertEquals(s2, start - 1)
        self.assertEquals(e2, s2+2)

    def test_pagging_token(self):
        request0 = self.factory.post(
            "/data/ga4gh/variants/search", json.dumps(
                {"end": 51425158029, "referenceName": "chr13", "variantSetId": "brca-hg37", "start": 51425158,
                 "pageSize": 5}), content_type="application/json")
        response = variant_search(request0)
        self.assertEqual(json.loads(response.content)["nextPageToken"], "1" )

        request0 = self.factory.post(
            "/data/ga4gh/variants/search", json.dumps(
                {"end": 51425158029, "referenceName": "chr17", "variantSetId": "brca-hg37", "start": 51425158,
                 "pageSize": 5, "pageToken": "3"}), content_type="application/json", )
        response = variant_search(request0)
        self.assertEqual(json.loads(response.content)["nextPageToken"], "4")
        self.assertEquals(len(json.loads(response.content)["variants"]), 5)

    def test_get_variant_response(self):
        request = self.factory.get("/data/ga4gh/variants/hg37-1")
        response = get_var_by_id(request, "hg37-1")
        jsonresp = json.loads(response.content)
        self.assertIsNotNone(jsonresp["referenceName"], True)
        self.assertIsNotNone(jsonresp["start"], True)

    def test_brca_to_ga4gh_variantSets_status_code(self):
        request0 = self.factory.post(
            "/data/ga4gh/variants/search", json.dumps({"datasetId": "hg37"}), content_type="application/json")
        response = get_variantSet(request0)
        self.assertEqual(response.status_code, 200)

    def test_brca_to_ga4gh_variantSets_metafieldNums(self):
        request0 = self.factory.post(
            "/data/ga4gh/variantsets/search", json.dumps({"datasetId": "brca-exchange"}), content_type="application/json")
        response = get_variantSet(request0)
        self.assertEqual(len(json.loads(response.content)["variantSets"][0]["metadata"]), 82)

    def test_variantSets_fields(self):
        request0 = self.factory.post(
            "/data/ga4gh/variants/search", json.dumps({"datasetId": "brca-exchange"}), content_type="application/json")
        response = get_variantSet(request0)
        Comp = json.loads(response.content)["variantSets"][1]

        """Comparative fields returned by request's responce"""
        self.assertEqual(Comp["referenceSetId"], "Genomic-Coordinate-hg37")
        self.assertEqual(Comp["datasetId"], "brca-exchange")
        self.assertEqual(Comp["name"], "brca-exchange-variants-hg37")

        """Next Page Token Field"""
        self.assertEquals(json.loads(response.content)["nextPageToken"], '0')


    def test_get_variantSet_by_id(self):
        request = self.factory.get("/data/ga4gh/variantsets/brca-hg37")
        response = get_varset_by_id(request, "brca-hg37")
        jsonresp = json.loads(response.content)
        self.assertIsNotNone(jsonresp["referenceSetId"], True)
        self.assertIsNotNone(jsonresp["id"], True)

    def test_empty_request(self):
        request0 = self.factory.get("data/ga4gh/variantsets")
        request1 = self.factory.get("data/ga4gh/variants")
        response0 = varsetId_empty_catcher(request0)
        response1 = empty_varId_catcher(request1)
        self.assertJSONEqual(response0.content,{'error code': 400, 'message' : 'invalid request empty request'})
        self.assertJSONEqual(response1.content, {'error code': 400, 'message': 'invalid request empty request'})

    def test_search_datasets(self):
        request = self.factory.post(
            "/data/ga4gh/datasets/search",  content_type="application/json")
        resoponse = search_datasets(request)
        jresponse = json.loads(resoponse.content)
        self.assertIsNotNone(jresponse["datasets"], True)
        self.assertIsNotNone(jresponse["datasets"][0], True)


if __name__ == '__main__':
    unittest.main()
