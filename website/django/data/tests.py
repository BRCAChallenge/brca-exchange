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
import data.views as views


# GA4GH related modules
import google.protobuf.json_format as json_format
from ga4gh import variant_service_pb2 as variant_service
from ga4gh import variants_pb2 as variants

from django.test import Client
c = Client()


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

    #######################
    # GA4GH related tests #
    #######################

    ######################
    # Datasets endpoints #
    ######################
    def test_search_datasets_alive(self):
        """Tests that datasets search endpoint is alive."""
        request = self.factory.post(
            "/data/ga4gh/datasets/search", json.dumps({"pageSize": 5}), content_type="application/json")
        response = views.search_datasets(request)
        self.assertEqual(response.status_code, 200)

    def test_search_datasets(self):
        """Ensures the responses from datasets search are well formed."""
        # TODO make sure the response is well formed here
        request = self.factory.post(
            "/data/ga4gh/datasets/search",  content_type="application/json")
        response = views.search_datasets(request)
        jresponse = json.loads(response.content)
        self.assertIsNotNone(jresponse["datasets"])
        self.assertIsNotNone(jresponse["datasets"][0])

    def test_search_datasets_paging(self):
        # TODO
        pass

    ##########################
    # Variant sets endpoints #
    ##########################
    def test_search_variant_sets_alive(self):
        """Tests that variantsets search endpoint is alive."""
        request = self.factory.post(
            "/data/ga4gh/variants/search",
            json.dumps({"datasetId": "brca"}),
            content_type="application/json")
        response = views.search_variant_sets(request)
        self.assertEqual(response.status_code, 200)

    def test_search_variant_sets(self):
        """Tests that variantsets search endpoint is alive."""
        # TODO show bad requests and good request
        pass

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

        # Next Page Token Field
        self.assertEquals(json.loads(response.content)["nextPageToken"], "0")

    def test_get_variant_set_by_id(self):
        """Ensures that expected variant sets are present."""
        # TODO get each variant set present here by ID
        request = self.factory.get("/data/ga4gh/variantsets/brca-hg37")
        response = views.get_variant_set(request, "brca-hg37")
        jsonresp = json.loads(response.content)
        self.assertIsNotNone(jsonresp["referenceSetId"], True)
        self.assertIsNotNone(jsonresp["id"], True)

    def test_variant_set_paging(self):
        # TODO page 1 by 1
        pass


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
        # TODO validate against mistyped entries, like if datasetid is not a string
        # TODO send a useful error if start > end
        request = variant_service.SearchVariantsRequest()
        req = self.factory.post("/data/ga4gh/variants/search",
                                json.dumps(json_format._MessageToJsonObject(request, False)),
                                content_type="application/json")
        response = views.search_variants(req)
        self.assertEqual(response.status_code, 400, "No variant set ID should 400")
        self.assertJSONEqual(response.content,
                             views.ErrorMessages['variantSetId'],
                             "No variant set ID in the request should provide a useful error")

        request.variant_set_id = "Something not null"
        req = self.factory.post("/data/ga4gh/variants/search",
                                json.dumps(json_format._MessageToJsonObject(request, False)),
                                content_type="application/json")
        response = views.search_variants(req)
        self.assertEquals(response.status_code, 400)
        self.assertJSONEqual(response.content,
                             views.ErrorMessages['referenceName'],
                             "A useful error is thrown when the reference name is not present")

        request.reference_name = "chr17"
        req= self.factory.post("/data/ga4gh/variants/search",
                               json.dumps(json_format._MessageToJsonObject(request, False)),
                                content_type="application/json")
        response = views.search_variants(req)
        self.assertJSONEqual(response.content,
                             views.ErrorMessages['start'],
                             "A useful error is thrown when no start is present")
        request.start = 14589
        req = self.factory.post("/data/ga4gh/variants/search",
                                json.dumps(json_format._MessageToJsonObject(request, False)),
                                content_type="application/json")
        response = views.search_variants(req)
        self.assertJSONEqual(response.content, views.ErrorMessages['end'],
                             "A useful error is provided when no end is present")

        request.end = 143295
        req = self.factory.post("/data/ga4gh/variants/search",
                                json.dumps(json_format._MessageToJsonObject(request, False)),
                                content_type="application/json")
        response = views.search_variants(req)
        self.assertEquals(response.status_code, 404, "A bad variant set ID should 404")

    def test_search_variants_paging_token(self):
        """Tests that paging works as expected for the search variants endpoint."""
        # Request a very large range of variants
        # TODO add a test that shows what happens when a known empty page is
        # requested
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
        start = 42245664
        end = 41245664
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
            self.assertGreater(variant["end"], start)
            self.assertLess(variant["start"], end,
                            "Should be in range"
                            " v.start {} r.end {}".format(variant['start'], end))
        self.assertEqual(response["nextPageToken"], "4")
        self.assertEquals(len(json.loads(response.content)["variants"]), 5)

    def test_search_variants_requested_range_present(self):
        """Ensures variants returned via search have the expected range."""
        # TODO add a test to show they're returned in order
        # TODO add a test to show when no variants are found in the range
        # the response is still well formed: {"variants": []}
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
            self.assertTrue(variant["end"] > start)
            self.assertTrue(variant["start"] < end,
                            "Should be in range"
                            " v.start {} r.end {}".format(variant['start'], end))
            self.assertTrue(variant["referenceName"] == "chr17")
        request = self.factory.post("/data/ga4gh/variants/search",
                                    json.dumps({"end": end,
                                                "referenceName": "17",
                                                "variantSetId": "brca-hg37",
                                                "start": start}),
                                    content_type = "application/json")
        response1 = views.search_variants(request)
        self.assertEqual(response.content, response1.content,
                         "Both 17 and chr17 return the same results")

    def test_get_variant_by_id(self):
        """Ensures the results found via search variants and get variant by ID are equal."""
        # TODO first search for a variant then get it
        # TODO test for a known bad ID gives 404
        request = self.factory.get("/data/ga4gh/variants/hg37-1")
        response = views.get_variant(request, "hg37-19")
        jsonresp = json.loads(response.content)
        self.assertIsNotNone(jsonresp["referenceName"])
        self.assertIsNotNone(jsonresp["start"])
        self.assertIsNotNone(jsonresp["end"])

    def test_empty_request(self):
        # TODO move to test search variants and sets respectively
        # TODO need same test for datasets
        request0 = self.factory.get("data/ga4gh/variantsets")
        request1 = self.factory.get("data/ga4gh/variants")
        response0 = views.varsetId_empty_catcher(request0)
        response1 = views.empty_varId_catcher(request1)
        self.assertJSONEqual(response0.content,
                             views.ErrorMessages['emptyBody'])
        self.assertJSONEqual(response1.content,
                             views.ErrorMessages['emptyBody'])

    def test_origin(self):
        # TODO tests to see that the HGVS and start strings represent the 1-based
        # and 0-based origins respectively
        pass

    def test_brca_to_ga4gh(self):
        """
        Gets rows from database and converts them into tests that
        ensure the function brca_to_ga4gh makes well formed GA4GH messages.
        """
        # TODO redesign this test
        test_variant = variants.Variant()
        db_variant = Variant.objects.values()[0]
        reference_genome = "hg37"
        for key in db_variant:
            if db_variant[key] == "-" or (db_variant[key] == ""):
                continue
            if key == "Genomic_Coordinate_" + reference_genome:
                # Parsing out position from genomic coordinate string
                bases = db_variant[key].split(':')[2]
                test_variant.reference_bases, alternbases = bases.split(">")
                # TODO is this correct?
                for i in range(len(alternbases)):
                    test_variant.alternate_bases.append(alternbases[i])
                continue
            if key == "id":
                test_variant.id = reference_genome + "-" + str(db_variant['id'])
                test_variant.variant_set_id = "brca" + "-" + reference_genome
                # TODO set a created an updated time based on something
                # This shouldn't pass like this
                test_variant.created = 0
                test_variant.updated = 0
                continue
            if key == "Reference_Name":
                test_variant.reference_name = db_variant[key]
                continue

            if key == "Hg36_Start" and (reference_genome == "hg36"):
                test_variant.start = db_variant[key]
                test_variant.end = db_variant["Hg36_End"]
                continue
            if key == "Hg36_End" and (reference_genome == "hg36"): continue
            if key == "Hg37_Start" and (reference_genome == "hg37"):
                test_variant.start = db_variant[key]
                test_variant.end = db_variant["Hg37_End"]
                continue
            if key == "Hg37_End" and (reference_genome == "hg37"): continue

            if key == "Hg38_Start" and (reference_genome == "hg38"):
                test_variant.start = db_variant[key]
                test_variant.end = db_variant["Hg38_End"]
                continue
            if key == "Hg38_End" and (reference_genome == "hg38"): continue

            if key == "Synonyms":
                names = [i for i in str(db_variant[key]).split(",")]
                for name in names:
                    test_variant.names.append(name)
            else:
                test_variant.info[str(key)].append(db_variant[key])

        expectedResp = json.dumps(json_format._MessageToJsonObject(test_variant, False))
        jresponse = json.dumps(
            json_format._MessageToJsonObject(
                views.brca_to_ga4gh(db_variant, reference_genome), False))
        self.assertJSONEqual(jresponse, expectedResp)

    def test_offset_calculator(self):
        """Tests that the ga4gh origin offset returns expected results."""
        start = 12500
        end = 13000
        # Ending offset is the same, the query will be exclusive, and
        # therefore query will satisfy type, that is inclusive.
        s1, e1 = views.get_offset(start, end)
        self.assertEquals(s1, start + 1)
        self.assertEquals(e1, end + 1)
        opperlen = ['a', 'b']
        s2, e2 = views.get_offset(start, end, opperlen)
        self.assertEquals(s2, start - 1)
        self.assertEquals(e2, s2 + 2)


if __name__ == '__main__':
    unittest.main()
