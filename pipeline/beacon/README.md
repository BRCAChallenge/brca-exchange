# Query Beacons for BRCA Variants

This script uses the GA4GH API to acquire a list of all variants present in the
BRCA Exchange, and then queries the Beacon Network for each. The output is a
JSON file containing data on which beacons have each variant, as well as data
on all the beacons currently in the network.

To run the script, run:

    query_beacons

Options can be viewed using:

    query_beacons --help

The GA4GH client package and the Requests package must be installed.
To install, run:

    pip install requests ga4gh==0.3.5
