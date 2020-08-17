'''
Parses original ClinVar XML file and filters out ClinVarSets of interest (i.e.
currently those related to BRCA1 and BRCA2)
'''

import gzip
import itertools

import click
from lxml import etree
from clinvar import clinvar_common

def open_maybe_gzip(fname):
    if fname.endswith('gz'):
        return gzip.GzipFile(fname)
    else:
        return open(fname, 'r')


def filter_xml(fin, fout, symbols):
    fout = open(fout, 'wb')

    xpath_filter = clinvar_common.build_xpath_filter_for_cv_assertions(symbols)

    # copying header
    with open_maybe_gzip(fin) as f:
        header_lines = itertools.takewhile(lambda s: "ClinVarSet" not in s.decode('utf-8'), f)
        fout.writelines(header_lines)

    # filtering ClinVarSet's we are interested in
    with open_maybe_gzip(fin) as f:
        for _, el in etree.iterparse(f, tag='ClinVarSet'):
            if len(el.xpath(xpath_filter)) >= 1:
                fout.write(etree.tostring(el, pretty_print=True, encoding='UTF-8'))

            # some cleanup to save a significant amount of memory
            # (inspired by https://www.ibm.com/developerworks/xml/library/x-hiperfparse/)
            el.clear()
            # Also eliminate now-empty references from the root node to elem
            for ancestor in el.xpath('ancestor-or-self::*'):
                while ancestor.getprevious() is not None:
                    del ancestor.getparent()[0]

    # writing "footer"
    fout.write("</ReleaseSet>".encode('UTF-8'))

    fout.close()


@click.command()
@click.argument('input', type=click.Path(exists=True))
@click.argument('output', type=click.Path())
@click.option('--gene', type=str, required=True, multiple=True)
def main(input, output, gene):
    filter_xml(input, output, list(set(gene)))

if __name__ == "__main__":
    main()
