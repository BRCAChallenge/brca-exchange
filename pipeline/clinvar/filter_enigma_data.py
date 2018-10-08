'''
Parses an Clinvar XML and filters out all ClinVarSet elements which contain ENIGMA submissions.
'''

import click
import lxml
from lxml import etree


def filter_engima_xml(fin, fout):
    f = etree.parse(fin)
    root = f.getroot()

    enigma_assertions = root.xpath(
        '//ClinVarSet/ClinVarAssertion/ClinVarSubmissionID[contains(@submitter, "ENIGMA")]')

    # careful: the more natural xpath expression seems (not needing getparent().getparent() later):
    #  //ClinVarSet[contains(ClinVarAssertion/ClinVarSubmissionID/@submitter, "ENIGMA")]
    # but this doesn't work as expected, as the enigma assertions are only picked up
    # if they are the first element.

    # copy ClinVarSet elements with enigma assertions to new document
    new_root = lxml.etree.Element(root.tag)
    [new_root.append(e.getparent().getparent()) for e in enigma_assertions];

    # writing out.
    etree.ElementTree(new_root).write_c14n(fout)


@click.command()
@click.argument('clinvar_xml', type=click.Path(exists=True))
@click.argument('output', type=click.Path())
def main(clinvar_xml, output):
    filter_engima_xml(clinvar_xml, output)
    # TODO: return code.


if __name__ == "__main__":
    main()
