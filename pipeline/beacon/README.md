# Query Beacons for BRCA Variants

This script uses the GA4GH API to acquire a list of all variants present in the BRCA Exchange, and then queries the Beacon Network for each. The output is a Python dictionary containing the Beacon Network website url for each variant, as well as a count of the number of beacons that responded affirmatively to the request.

To run the script, run:

    python query_beacons.py

The GA4GH client package and the Requests package must be installed. To install, run:

    pip install requests ga4gh-client
