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
        return open(fname)


def filter_engima_xml(fin, fout):
    fout = open(fout, 'w')

    xpath_filter = 'ClinVarSet/ReferenceClinVarAssertion/MeasureSet/Measure/' \
                   'MeasureRelationship/Symbol[starts-with(ElementValue, "BRCA") and ElementValue/@Type="Preferred"]'

    # writing header
    with open_maybe_gzip(fin) as f:
        header_lines = itertools.takewhile(lambda s: "ClinVarSet" not in s, f)
        fout.writelines(header_lines)

    # filtering ClinVarSet's we are interested in
    for ev_name, el in etree.iterparse(fin, tag='ClinVarSet'):
        if not el.xpath(xpath_filter):
            fout.write(etree.tostring(el))

    # writing "footer"
    fout.write("</ReleaseSet>")

    fout.close()


@click.command()
@click.argument('input', type=click.Path(exists=True))
@click.argument('output', type=click.Path())
def main(input, output):
    filter_engima_xml(input, output)


if __name__ == "__main__":
    main()
