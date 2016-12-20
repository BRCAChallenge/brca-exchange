"""
Author:Audrey Musselman-Brown
"""

import sys
import requests
import json
import time
import threading
from ga4gh.client import client

URL_MAX = 1000  # the longest allowed URL length
TOTAL_TRIES = 1 # the number of times to attempt each beacon request
SLEEP_TIME = 1  # the number of seconds to sleep between requests

VariantD = {}

beacon_api_url = "https://beacon-network.org/api/responses/"

try:
    beacons = json.loads(requests.get("https://beacon-network.org/api/beacons").content)
except:
    sys.exit("Failed to access the Beacon Network")


class Gene (object):

    def __init__(self, chromosome, start, end):
        object.__init__(self)
        self.chromosome = chromosome
        self.start = start
        self.end = end

brca1 = Gene("chr17", 41160094-5 , 41322380+5)
brca2 = Gene("chr13", 32889080-5, 32973779+5)

class QueryBeacons (threading.Thread):
    """
    A thread object that requests the given variant from all beacons
    """
    def __init__(self, variant):
        threading.Thread.__init__(self)
        self.variant = variant.id
        self.request = "?chrom={}&pos={}&allele={}&ref=GRCh37".format(
                          variant.reference_name, variant.start,
                          variant.alternate_bases[0])

    def run(self):

        for beacon in beacons:
            
            if beacon["id"] == "bob":
                continue
            
            success = False
            tries = 0
            while not success and tries <= 5:
                try:
                    r = requests.get(beacon_api_url + beacon["id"] + self.request)
                    content = json.loads(r.content)
                    if "response" in content:
                        if content["response"] is not None:
                            VariantD[self.variant]["beacons with variant"] += bool(content["response"])
                            VariantD[self.variant]["beacon responses"] += 1
                    else:
                        sys.stderr.write("no response from "+beacon["id"]+"\n")
                    success = True

                except:
                    tries += 1
                    if tries >= TOTAL_TRIES:
                        print "Failed request: "+beacon["id"]
                        raise


def main(args):

    if len(args) > 0 and args[0] == "--test":
        brca1.end = brca1.start + 100
        brca2.end = brca2.start + 100

    c = client.HttpClient("http://brcaexchange-dev.cloudapp.net:9004/data/ga4gh/v0.6.0a7/")

    threads = []

    variants_so_far = 0
    for gene in [brca1, brca2]:
        for variant in c.search_variants(reference_name=gene.chromosome,
                                         variant_set_id="brca-hg37",
                                         start=gene.start,
                                         end=gene.end):
            #print variant.id
            #print variant.variant_set_id
            #print variant.reference_name
            #print variant.start
            #print variant.end
            #print variant.reference_bases
            #print variant.alternate_bases
            
	    in_enigma = variant.info[u'Variant_in_ENIGMA'][0]
            
            if in_enigma:
                pass

            if variants_so_far % 100 == 0:
                sys.stderr.write(str(i) + " variants processed\n")
            variants_so_far += 1
           
            website_url = "https://beacon-network.org/#/search?chrom={}&pos={}&ref={}&allele={}&rs=GRCh37".format(
                              variant.reference_name, variant.start,
                              variant.reference_bases, variant.alternate_bases[0])

            if len(request) > URL_MAX:
                sys.stderr.write("variant too long\n")
                continue

            VariantD[variant.id] = {"url": website_url, "beacons with variant": 0, "beacon responses": 0}

            thread = QueryBeacons(variant)
            thread.start()
            threads.append(thread)
            time.sleep(SLEEP_TIME)
        
    for thread in threads:

        thread.join()

    print VariantD

if __name__=="__main__":
    sys.exit(main(sys.argv[1:]))

