#!/usr/bin/env python
import datetime
import json
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--date", default=datetime.date.today(),
                        help="Version generation date")
    parser.add_argument("--notes", required=True,
                        help="File with release notes text")
    parser.add_argument("--output", default="version.json", help="Output json file")

    args = parser.parse_args()
    date = args.date

    version_data = {}
    version_data['date'] = date

    # Make sure this matches the format of the release archive created from the pipeline in GenerateReleaseArchive.
    date = datetime.datetime.strptime(date, "%Y-%m-%d")
    version_data['archive'] = "release-" + date.strftime("%x").replace('/', '-') + ".tar.gz"

    release_notes = None
    with open(args.notes) as release_notes_file:
        release_notes = release_notes_file.read()

    version_data['notes'] = release_notes

    # TODO: come up with a way to generate this and creation dates programatically.
    source_files = ["Bic", "ClinVar", "ESP", "ExAC", "ENIGMA", "LOVD", "ExUV", "1000 Genomes"]
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
