import json
import os
import unittest
from unittest import skip
from urllib import quote
from django.http import JsonResponse, HttpResponse
from django.test import TestCase
from django.test.client import RequestFactory
from brca import settings
from data import test_data
from data.models import Variant
from data.views import index, autocomplete
import data.views as views

# GA4GH related modules
import google.protobuf.json_format as json_format
from ga4gh.schemas.ga4gh import variant_service_pb2 as variant_service

class VariantTestCase(TestCase):
    def setUp(self):
        self.factory = RequestFactory()
        datafile = os.path.join(settings.BASE_DIR, 'data', 'resources', 'releases', 'release-10-06-16', 'built_with_change_types.tsv')
        self.db_size = sum(1 for _ in open(datafile)) - 1

    @skip("not complete")
    def test_variant_model(self):
        """Create a new variant and then retrieve it by the Genomic_Coordinate_hg38 column"""
        self.assertEqual(len(Variant.objects.all()), self.db_size)
        Variant.objects.create_variant(row=(test_data.new_variant()))
        self.assertEqual(len(Variant.objects.all()), self.db_size + 1)
        retrieved_variant = Variant.objects.get(Genomic_Coordinate_hg38="chr17:999999:A>G")
        self.assertIsNotNone(retrieved_variant)

    @skip("Not Complete")
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

    @skip("Not complete")
    def test_autocomplete_nucleotide(self):
        """Getting autocomplete suggestions for words starting with c.2123 should return 2 results"""
        search_term = quote('c.2123')
        expected_autocomplete_results = [["c.2123c>a"], ["c.2123c>t"]]

        request = self.factory.get('/data/suggestions/?term=%s' % search_term)
        response = autocomplete(request)

        self.assertIsInstance(response, JsonResponse)
        self.assertEqual(response.status_code, 200)
        self.assertJSONEqual(response.content, {"suggestions": expected_autocomplete_results})

    @skip("Not complete")
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

    #######################
    # GA4GH related tests #
    #######################

    ######################
    # Datasets endpoints #
    ######################
    def test_datasets_alive(self):
        """Tests that datasets search endpoint is alive."""
        request = self.factory.post(
            "/data/ga4gh/datasets/search",
            json.dumps({"pageSize": 1}),
            content_type="application/json")
        response = views.search_datasets(request)
        self.assertEqual(response.status_code, 200)
        request_2 = self.factory.get(
            "/data/ga4gh/datasets/brca")
        response_2 = views.get_dataset(request_2, "brca")
        self.assertEqual(response_2.status_code, 200)

    def test_search_datasets(self):
        """Ensures the responses from datasets search are well formed."""
        request = self.factory.post(
            "/data/ga4gh/datasets/search",
            content_type="application/json")
        response = views.search_datasets(request)
        json_response = json.loads(response.content)
        self.assertIsNotNone(json_response["datasets"])
        for dataset in json_response["datasets"]:
            self.assertIsNotNone(dataset)

    def test_search_datasets_paging(self):
        """Server dataset field should be empty, since server only supports one dataset """
        request = self.factory.post(
            "/data/ga4gh/datasets/search",
            json.dumps({"pageSize": 1, "nextPageToken": "2"}),
            content_type="application/json")
        response = views.search_datasets(request)
        json_response = json.loads(response.content)
        # Note that it was asserted to an empty string, because it is not a none value.
        self.assertTrue(json_response.get("nextPageToken", None) is None)

    def test_get_datasets(self):
        request_2 = self.factory.get(
            "/data/ga4gh/datasets/brca")
        response_2 = views.get_dataset(request_2, "brca")
        json_response = json.loads(response_2.content)
        self.assertIsNotNone(json_response)
        for elements in json_response:
            self.assertIsNotNone(elements)

    def test_empty_datasets_request(self):
        request = self.factory.get("data/ga4gh/datasets")
        response = views.empty_variantset_id_catcher(request)
        self.assertEqual(response.status_code, 405)
        self.assertJSONEqual(response.content,
                             views.ErrorMessages['methodNotAllowed'])

    ##########################
    # Variant sets endpoints #
    ##########################
    def test_search_and_get_variant_sets_alive(self):
        """Tests that variantsets search endpoint is alive."""
        search_request = self.factory.post(
            "/data/ga4gh/variantsets/search",
            json.dumps({"datasetId": "brca"}),
            content_type="application/json")
        response = views.search_variant_sets(search_request)
        self.assertEqual(response.status_code, 200)
        get_request_hg36 = self.factory.get(
            "/data/ga4gh/variantsets/brca-hg36")
        response_hg36 = views.get_variant_set(get_request_hg36, "brca-hg36")
        self.assertEqual(response_hg36.status_code, 200)

        get_request_hg37 = self.factory.get(
            "/data/ga4gh/variantsets/brca-hg37")
        response_hg37 = views.get_variant_set(
            get_request_hg37, "brca-hg37")
        self.assertEqual(response_hg37.status_code, 200)

        get_request_hg38 = self.factory.get(
            "/data/ga4gh/variantsets/brca-hg38")
        response_hg38 = views.get_variant_set(
            get_request_hg38, "brca-hg38")
        self.assertEqual(response_hg38.status_code, 200)

    def test_search_variant_sets(self):
        """Tests that variantsets search endpoint is alive."""
        bad_set_id_parameter = "not-brca1002"
        search_request = self.factory.post(
            "/data/ga4gh/variantsets/search",
            json.dumps({"datasetId": bad_set_id_parameter}),
            content_type="application/json")
        response = views.search_variant_sets(search_request)
        self.assertEqual(response.status_code, 200)
        # Bad request 'data_set_id' returns a valid, empty response
        json_response = json.loads(response.content)
        """Note that the following test are assertEquals do to the response,
        not being empty, therefore we check for empty values."""
        # Field 'callSets' is an empty list in bad request set id
        self.assertEqual(json_response["callSets"], list([]))
        # Next page token is empty, no other result available for display
        self.assertEqual(json_response["nextPageToken"], unicode(''))

    def test_variant_set_metadata(self):
        """Tests that a variant set found via search has well formed metadata."""
        request = self.factory.post(
            "/data/ga4gh/variantsets/search",
            json.dumps({"datasetId": "brca"}), content_type="application/json")
        response = views.search_variant_sets(request)
        metadata = json.loads(response.content)["variantSets"][0]["metadata"]
        for item in metadata:
            self.assertIsNotNone(item.get('id', None))

    def test_variant_set_fields(self):
        """Tests that variant set messages returned via search are well formed."""
        request = self.factory.post(
            "/data/ga4gh/variants/search",
            json.dumps({"datasetId": "brca"}),
            content_type="application/json")
        response = views.search_variant_sets(request)
        variant_sets = json.loads(response.content)["variantSets"]
        # There are currently three variant sets
        self.assertEquals(len(variant_sets), 3)
        for variant_set in variant_sets:
            self.assertIsNotNone(variant_set["referenceSetId"])
            self.assertIsNotNone(variant_set["datasetId"])
            self.assertIsNotNone(variant_set["name"])
        # Next Page Token Field should be empty response.
        self.assertEquals(json.loads(response.content)["nextPageToken"], "")

    def test_get_variant_set_by_id(self):
        """Ensures that expected variant sets are present."""
        # Test for variant set brca-hg36
        request_hg36 = self.factory.get("/data/ga4gh/variantsets/brca-hg36")
        response_hg36 = views.get_variant_set(request_hg36, "brca-hg36")
        json_response36 = json.loads(response_hg36.content)
        self.assertIsNotNone(json_response36["referenceSetId"])
        self.assertIsNotNone(json_response36["id"])
        # Test for variant set brca-hg37
        request_hg37 = self.factory.get("/data/ga4gh/variantsets/brca-hg37")
        response_hg37 = views.get_variant_set(request_hg37, "brca-hg37")
        json_response37 = json.loads(response_hg37.content)
        self.assertIsNotNone(json_response37["referenceSetId"])
        self.assertIsNotNone(json_response37["id"])
        # Test for variant set brca-hg38
        request_hg38 = self.factory.get("/data/ga4gh/variantsets/brca-hg38")
        response_hg38 = views.get_variant_set(request_hg38, "brca-hg38")
        json_response38 = json.loads(response_hg38.content)
        self.assertIsNotNone(json_response38["referenceSetId"])
        self.assertIsNotNone(json_response38["id"])

    def test_variant_set_paging(self):
        search_request = self.factory.post(
            "/data/ga4gh/variantsets/search",
            json.dumps({"datasetId": "brca",
                        "pageSize": 1,
                        "pageToken": "0"}),
            content_type="application/json")
        response = views.search_variant_sets(search_request)
        json_response = json.loads(response.content)
        response_variants = json_response["variantSets"][0]
        self.assertEqual(json_response["nextPageToken"], "1")
        self.assertEqual(response_variants["id"], "brca-hg36")
        self.assertEqual(response_variants["datasetId"], "brca")
        self.assertEqual(response_variants["referenceSetId"], "Genomic-Coordinate-hg36")

        search_request = self.factory.post(
            "/data/ga4gh/variantsets/search",
            json.dumps({"datasetId": "brca",
                        "pageSize": 1,
                        "pageToken": "1"}),
            content_type="application/json")
        response = views.search_variant_sets(search_request)
        json_response = json.loads(response.content)
        response_variants = json_response["variantSets"][0]
        self.assertEqual(json_response["nextPageToken"], "2")
        self.assertEqual(response_variants["id"], "brca-hg37")
        self.assertEqual(response_variants["datasetId"], "brca")
        self.assertEqual(response_variants["referenceSetId"], "Genomic-Coordinate-hg37")

        search_request = self.factory.post(
            "/data/ga4gh/variantsets/search",
            json.dumps({"datasetId": "brca",
                        "pageSize": 1,
                        "pageToken": "2"}),
            content_type="application/json")
        response = views.search_variant_sets(search_request)
        json_response = json.loads(response.content)
        response_variants = json_response["variantSets"][0]
        self.assertEqual(json_response["nextPageToken"], "")
        self.assertEqual(response_variants["id"], "brca-hg38")
        self.assertEqual(response_variants["datasetId"], "brca")
        self.assertEqual(response_variants["referenceSetId"], "Genomic-Coordinate-hg38")

    def test_empty_variant_sets_request(self):
        request = self.factory.get("data/ga4gh/variantsets")
        response = views.empty_variantset_id_catcher(request)
        self.assertEqual(response.status_code, 405)
        self.assertJSONEqual(response.content,
                             views.ErrorMessages['methodNotAllowed'])
    ######################
    # Variants endpoints #
    ######################
    def test_search_variants_alive(self):
        """Ensures that search variants responds successfully for a known search"""
        request =  self.factory.post(
            "/data/ga4gh/variants/search",
            json.dumps(
                {"referenceName": "chr17",
                 "variantSetId": "brca-hg37",
                 "start": 51425158,
                 "end": 515158029,
                 "pageSize": 5 }), content_type="application/json")
        response = views.search_variants(request)
        self.assertEqual(response.status_code, 200)

    def test_search_variants_request_validation(self):
        """Ensures the search variants endpoint responds with expected failure modes"""
        request = variant_service.SearchVariantsRequest()
        req = self.factory.post("/data/ga4gh/variants/search",
                                json.dumps(json_format.MessageToDict(request, False)),
                                content_type="application/json")
        response = views.search_variants(req)
        self.assertEqual(response.status_code, 400, "No variant set ID should 400")
        self.assertJSONEqual(response.content,
                             views.ErrorMessages['variantSetId'],
                             "No variant set ID in the request should provide a useful error")
        request.variant_set_id = "Something not null"
        req = self.factory.post("/data/ga4gh/variants/search",
                                json.dumps(json_format.MessageToDict(request, False)),
                                content_type="application/json")
        response = views.search_variants(req)
        self.assertEquals(response.status_code, 400)
        self.assertJSONEqual(response.content,
                             views.ErrorMessages['referenceName'],
                             "A useful error is thrown when the reference name is not present")
        request.reference_name = "chr17"
        req= self.factory.post("/data/ga4gh/variants/search",
                               json.dumps(json_format.MessageToDict(request, False)),
                                content_type="application/json")
        response = views.search_variants(req)
        self.assertJSONEqual(response.content,
                             views.ErrorMessages['start'],
                             "A useful error is thrown when no start is present")
        request.start = 14589
        req = self.factory.post("/data/ga4gh/variants/search",
                                json.dumps(json_format.MessageToDict(request, False)),
                                content_type="application/json")
        response = views.search_variants(req)
        self.assertJSONEqual(response.content, views.ErrorMessages['end'],
                             "A useful error is provided when no end is present")
        request.end = 143295
        req = self.factory.post("/data/ga4gh/variants/search",
                                json.dumps(json_format.MessageToDict(request, False)),
                                content_type="application/json")
        response = views.search_variants(req)
        self.assertEquals(response.status_code, 404, "A bad variant set ID should 404")
        # Test for an end value less than the end value

        test_request = self.factory.post("/data/ga4gh/variants/search",
                                json.dumps({"referenceName": "chr17", "variantSetId": "brca-hg37", "start": 10000, "end": 1000}),
                                content_type="application/json")
        response_x = views.search_variants(test_request)
        self.assertEqual(response_x.status_code, 400)
        self.assertJSONEqual(response_x.content, views.ErrorMessages['invalidPositions'])

    @skip("Not complete")
    def test_search_variants_paging_token(self):
        """Tests that paging works as expected for the search variants endpoint."""
        # Request a very large range of variants
        request = self.factory.post(
            "/data/ga4gh/variants/search", json.dumps(
                {"end": 45158029,
                 "referenceName": "chr13",
                 "variantSetId": "brca-hg37",
                 "start": 20425158,
                 "pageSize": 5}), content_type="application/json")
        response = views.search_variants(request)
        self.assertEqual(
            json.loads(response.content)["nextPageToken"],
            "1", "A page token should be set")
        start = 41245664
        end = 42245664
        request = self.factory.post(
            "/data/ga4gh/variants/search", json.dumps(
                {"end": end,
                 "referenceName": "chr17",
                 "variantSetId": "brca-hg37",
                 "start": start,
                 "pageSize": 5,
                 "pageToken": "3"}), content_type="application/json", )
        response = json.loads(views.search_variants(request).content)
        for variant in response["variants"]:
            self.assertGreater(long(variant["end"]), start)
            self.assertLess(long(variant["start"]), end,
                            "Should be in range"
                            " v.start {} r.end {}".format(variant['start'], end))
        self.assertEqual(response["nextPageToken"], "4")
        self.assertEquals(len(response["variants"]), 5)

    def test_search_variants_empty_page_request(self):
        """This would be the same response when a requested
        range is not supported, that is, no variants are registered
        within the requested range"""

        request = self.factory.post(
            "/data/ga4gh/variants/search",
            json.dumps({"variantSetId": "brca-hg37",
                        "referenceName": "chr17",
                        "start": 10000, "end": 41196822,
                        "pageSize": 10, "pageToken": "5"}),
            content_type="application/json")
        response = views.search_variants(request)
        json_response = json.loads(response.content)
        # Note that it is an empty response, therefore assertions are done to empty values
        self.assertEqual(json_response["nextPageToken"], "")
        self.assertEqual(json_response["variants"], list([]))

    @skip("Not complete")
    def test_search_variants_requested_range_present(self):
        """Ensures variants returned via search have the expected range."""
        start = 41246794
        end = 41247374
        request = self.factory.post("/data/ga4gh/variants/search",
                                    json.dumps({"end": end,
                                                "referenceName": "chr17",
                                                "variantSetId": "brca-hg37",
                                                "start": start}),
                                    content_type = "application/json")
        response = views.search_variants(request)
        jresp = json.loads(response.content)
        self.assertGreater(len(jresp["variants"]), 0)
        for variant in jresp["variants"]:
            self.assertTrue(long(variant["end"]) > start)
            self.assertTrue(long(variant["start"]) < end,
                            "Should be in range"
                            " v.start {} r.end {}".format(variant['start'], end))
            self.assertTrue(
                variant["referenceName"] == "chr17" or variant["referenceName"] == "17",
                    "Searched for chr17 and got {}".format(variant["referenceName"]))
        request = self.factory.post("/data/ga4gh/variants/search",
                                    json.dumps({"end": end,
                                                "referenceName": "17",
                                                "variantSetId": "brca-hg37",
                                                "start": start}),
                                    content_type = "application/json")
        response1 = views.search_variants(request)
        self.assertEqual(response.content, response1.content,
                         "Both 17 and chr17 return the same results")

    def test_ordered_range_for_search_variant(self):
        START = 41246925
        END = 41247084
        # Assure Start is less than End
        self.assertGreater(END, START)
        request = self.factory.post("/data/ga4gh/variants/search",
                                    json.dumps({"referenceName": "chr17",
                                                "variantSetId": "brca-hg37",
                                                "start": START, "end": END,
                                                "pageSize": 9}),
                                    content_type="application/json")
        response = views.search_variants(request)
        json_response = json.loads(response.content)["variants"]
        previews_start = 0
        for variant in json_response:
            if previews_start == 0:
                previews_start = variant['start']
                self.assertGreaterEqual(variant["start"], START)
            else:
                """Note that if in order ranged is returned
                the previews start value should be less than the following
                in the list"""
                if variant['end'] == END:
                    break
                self.assertLess(previews_start, variant['start'])
                previews_start = variant['start']

    def test_search_unsupported_region(self):
        request = self.factory.post("/data/ga4gh/variants/search",
                                           json.dumps({"referenceName": "chr13",
                                                       "variantSetId": "brca-hg36",
                                                       "start": 10, "end": 10000}),
                                           content_type="application/json")
        response = views.search_variants(request)
        json_response = json.loads(response.content)
        """Note that it is an empty response, no variants exist within
        the provided search range"""
        self.assertEqual(json_response['nextPageToken'], unicode(''))
        self.assertEqual(json_response['variants'], list([]))
        request = self.factory.post("/data/ga4gh/variants/search",
                                           json.dumps({"referenceName": "5",
                                                       "variantSetId": "brca-hg36",
                                                       "start": 10, "end": 10000000000}),
                                           content_type="application/json")
        response = views.search_variants(request)
        json_response = json.loads(response.content)
        """Note that it is an empty response, no variants exist within
        the provided search range"""
        self.assertEqual(response.status_code, 200)
        self.assertEqual(json_response['nextPageToken'], unicode(''))
        self.assertEqual(json_response['variants'], list([]))

    @skip("Not complete")
    def test_get_variant_by_id(self):
        """Ensures the results found via search variants and get variant by ID are equal."""
        search_request = self.factory.post("/data/ga4gh/variants/search",
                                           json.dumps({"referenceName": "chr17",
                                                       "variantSetId": "brca-hg37",
                                                       "start" : 4124692, "end" : 41247086,
                                                       "pageSize": 1}),
                                           content_type="application/json"
                                           )
        search_response = views.search_variants(search_request)
        json_search_response = json.loads(search_response.content)["variants"][0]
        self.assertIsNotNone(json_search_response["variantSetId"])
        self.assertIsNotNone(json_search_response["referenceName"])
        self.assertIsNotNone(json_search_response["id"])

        search_based_id = str(json_search_response['id'])
        get_request = self.factory.get("/data/ga4gh/variants/"+search_based_id)
        response = views.get_variant(get_request, search_based_id)
        json_get_response = json.loads(response.content)

        """Note, because we made a get_request based on the search above,
        we should be able to confirm that the individual variant request is exactly the same one"""
        self.assertEqual(json_get_response["referenceName"], json_search_response["referenceName"])
        self.assertEqual(json_get_response["start"], json_search_response["start"])
        self.assertEqual(json_get_response["end"], json_search_response["end"])

    def test_bad_get_variant_request(self):
        invalid_id = "Jijf30453hwbur-PWFvWIvwPRGrjnib"
        request = self.factory.get("/data/ga4gh/variants/"+invalid_id)
        response = views.get_variant(request, invalid_id)
        # Useful error message for an id which is not supported
        self.assertEqual(response.status_code, 404)
        self.assertJSONEqual(response.content, views.ErrorMessages['notFoundId'])

        invalid_id_2 = "hg37-999999999999"
        request = self.factory.get("/data/ga4gh/variants/" + invalid_id_2)
        response = views.get_variant(request, invalid_id_2)
        # Useful error for elements non existent in data base"""
        self.assertEqual(response.status_code, 404)
        self.assertJSONEqual(response.content, views.ErrorMessages['notFoundId'])

    def test_empty_variants_request(self):
        request = self.factory.get("data/ga4gh/variants")
        response = views.empty_variant_id_catcher(request)
        self.assertEqual(response.status_code, 405)
        self.assertJSONEqual(response.content,
                             views.ErrorMessages['methodNotAllowed'])

    @skip("Not complete")
    def test_origin(self):
        variant_id = "hg37-55"
        request = self.factory.get("/data/ga4gh/variants/"+variant_id)
        response = views.get_variant(request, variant_id)
        json_response = json.loads(response.content)
        hgvs_start = int(
            json_response['info']["Genomic_Coordinate_hg37"][0].split(':')[1].replace('g.',''))
        # Note that `Genomic Coordinate hg37` is 1 base offset from its determined 'start' position"""
        self.assertEqual(int(json_response['start']), hgvs_start,
                         "`base_1_start`, obtained from 'Genomic Coordinate hg37' minus 1 should be equal to the `start` Variant parameter {} {}".format(
                             int(json_response['start']), hgvs_start - 1))

        variant_id = "hg38-945"
        request = self.factory.get("/data/ga4gh/variants/" + variant_id)
        response = views.get_variant(request, variant_id)
        json_response = json.loads(response.content)
        # Note that `start` is 0 base offset from its determined 'Genomic Coordinate' field"""
        hgvs_start = int(json_response['info']['Genomic_Coordinate_hg38'][0].split(':')[1].replace('g.',''))
        self.assertEqual(int(json_response['start']), hgvs_start - 1)

    def test_brca_to_ga4gh(self):
        """
        Gets rows from database and converts them into tests that
        ensure the function brca_to_ga4gh makes well formed GA4GH messages.
        """
        variant = {'Genomic_Coordinate_hg37': 'chr13:32923951:CA>C', 'Chr': '13',
                              'id': 1, 'Hg37_End': 32923951, 'Genomic_Coordinate_hg36': 'chr13:31821951:CA>C',
                              'Hg36_Start': 31821950, 'Hg37_Start': 32923950, 'Genomic_Coordinate_hg38': 'chr13:32349814:CA>C',
                              'Hg38_End': 32349814,'Hg36_End': 31821951, 'Hg38_Start': 32349813, 'Synonyms': 'U43746.1:c.7234+2936delA'}

        genomic_coordinate = "hg37"

        response = views.brca_to_ga4gh(variant, genomic_coordinate)

        json_response = json_format.MessageToDict(response, True)
        self.assertEqual(int(json_response['start']), 32923950)
        self.assertEqual(json_response['referenceBases'], "CA")
        self.assertEqual(json_response['alternateBases'][0], "C")
        self.assertEqual(json_response['referenceName'], "13")
        self.assertEqual(json_response['id'], "hg37-1")

if __name__ == '__main__':
    unittest.main()
