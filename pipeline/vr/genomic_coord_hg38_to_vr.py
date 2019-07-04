#!/usr/bin/env python
"""
Translate BRCA Exchange genomic variant IDs into GA4GH VR objects.

This assumes the availabiity or a SeqRepo REST API, running at the default
URL or elsewhere.
"""
import argparse
import csv
import io
import json
import VR
    
def main():
    args = parse_args()
    vr_mapper = VR.VRUtils(args.seqrepo_rest_service_url)
    output_file = io.open(args.output_file, 'w')
    output_file.write("Variant_id\tVR\n")
    input_file = csv.DictReader(open(args.input_file), delimiter="\t")
    for row in input_file:
        variant_id = row[args.key_column]
        variant_cDNA = row[args.pyhgvs_cDNA_column]
        if args.verbose:
            print("working on", variant_id, variant_cDNA)
        vr_variant = vr_mapper.to_vr(variant_cDNA, args.verbose)
        if vr_variant is not None:
            if args.verbose:
                print("translated to", json.dumps(vr_variant.as_dict()))
            output_file.write("%s\t%s\n" % (variant_id,
                                            json.dumps(vr_variant.as_dict())))
    output_file.close()
        
        

def parse_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file",
                        help="File containing the list of variants")
    parser.add_argument("-k", "--key_column",
                        default="pyhgvs_Genomic_Coordinate_38",
                        help="column in the input file with the variant ID")
    parser.add_argument("-c", "--pyhgvs_cDNA_column",
                        default="pyhgvs_cDNA",
                        help="Preferred HGVS cDNA column")
    parser.add_argument("-o", "--output_file",
                        help="File mapping the variant IDs to VR objects")
    parser.add_argument("-s", "--seqrepo_rest_service_url",
                        default="http://localhost:5000/seqrepo",
                        help="URL for the SeqRepo REST API")
    parser.add_argument("-v", "--verbose", default=False,
                        help="Whether or not to show verbose output")
    args = parser.parse_args()
    return(args)

    
if __name__ == "__main__":
    main()
