'''
Parses original ClinVar XML file and filters out ClinVarSets of interest (i.e.
currently those related to BRCA1 and BRCA2)
'''

import gzip
import itertools

import click
from lxml import etree


def open_maybe_gzip(fname):
    if fname.endswith('gz'):
        return gzip.open(fname, 'rb')
    else:
        return gzip.GzipFile(fname)


def filter_xml(fin, fout):
    fout = open(fout, 'w')

    # copying ClinVarSet if it contains a Symbol with a "Preferred" ElementValue BRCA*
    xpath_filter = 'ReferenceClinVarAssertion/MeasureSet/Measure/' \
                   'MeasureRelationship/Symbol[starts-with(ElementValue, "BRCA") and ElementValue/@Type="Preferred"]'

    # copying header
    with open_maybe_gzip(fin) as f:
        header_lines = itertools.takewhile(lambda s: "ClinVarSet" not in s, f)
        fout.writelines(header_lines)

    # filtering ClinVarSet's we are interested in
    with open_maybe_gzip(fin) as f:
        for _, el in etree.iterparse(f, tag='ClinVarSet'):
            if len(el.xpath(xpath_filter)) >= 1:
                fout.write(etree.tostring(el, pretty_print=True))

            # some cleanup to save a significant amount of memory
            # (inspired by https://www.ibm.com/developerworks/xml/library/x-hiperfparse/)
            el.clear()
            # Also eliminate now-empty references from the root node to elem
            for ancestor in el.xpath('ancestor-or-self::*'):
                while ancestor.getprevious() is not None:
                    del ancestor.getparent()[0]

    # writing "footer"
    fout.write("</ReleaseSet>")

    fout.close()


@click.command()
@click.argument('input', type=click.Path(exists=True))
@click.argument('output', type=click.Path())
def main(input, output):
    filter_xml(input, output)


if __name__ == "__main__":
    main()
