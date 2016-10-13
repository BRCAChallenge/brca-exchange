#!/usr/bin/env python
import datetime
import json
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--date", default=datetime.datetime.now(),
                        help="Version generation date")
    parser.add_argument("--notes", required=True,
                        help="File with release notes text")
    parser.add_argument("--output", default="version.json", help="Output json file")

    args = parser.parse_args()

    version_data = {}
    version_data['date'] = args.date

    release_notes = None
    with open(args.notes) as release_notes_file:
        release_notes = release_notes_file.read()

    version_data['notes'] = release_notes

    # TODO: come up with a way to generate this and creation dates programatically.
    source_files = ["Bic", "ClinVar", "ESP", "ExAC", "Enigma", "LOVD", "exLOVD", "1000 Genomes"]
    version_data['sources'] = source_files

    json_data = json.dumps(version_data, default=handler)

    with open(args.output, 'w') as json_output:
        json_output.write(json_data)


def handler(obj):
    if hasattr(obj, 'isoformat'):
        return obj.isoformat()
    else:
        raise TypeError, 'Object of type %s with value of %s is not JSON serializable' % (type(obj), repr(obj))

if __name__ == "__main__":
    main()
