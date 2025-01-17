'''
Parses an Clinvar XML and filters out all ClinVarSet elements which contain ENIGMA submissions.
'''

import click
import re
import xml.etree.ElementTree as ET

def filter_engima_xml(fin, fout):
    context = ET.iterparse(fin, events=("start", "end"))
    (event, root) = next(context)
    new_root = ET.Element(root.tag, attrib=root.attrib)
    new_tree = ET.ElementTree(new_root)
    for (event, elem) in context:
        if event == 'end' and elem.tag == "VariationArchive":
            enigma_variant = False
            for sub_elem in elem.findall(".//ClinicalAssertion/ClinVarAccession"):
                if "SubmitterName" in sub_elem.attrib:
                    if re.search("ENIGMA", sub_elem.attrib["SubmitterName"]):
                        enigma_variant = True
            if enigma_variant:
                new_root.append(elem)
            else:
                elem.clear()
    # writing out.
    ET.ElementTree(new_root).write(fout)
    
@click.command()
@click.argument('clinvar_xml', type=click.Path(exists=True))
@click.argument('output', type=click.Path())
def main(clinvar_xml, output):
    filter_engima_xml(clinvar_xml, output)


if __name__ == "__main__":
    main()
