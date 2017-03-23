#!/usr/bin/env python
"""
Given an ENIGMA ClinVar spreadsheet, extract the contents in TSV format, minus the extra headers,
and with a couple key changes in the labels
"""
import argparse
import openpyxl
import re
import sys
reload(sys)
sys.setdefaultencoding('utf8')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_excel_file",
                        help="Excel file containing the submitted variant data")
    parser.add_argument("-o", "--output_tsv_file",
                        help="Output file in TSV format")
    args = parser.parse_args()

    tsv = open(args.output_tsv_file, "wb")
    xlsx = openpyxl.load_workbook(args.input_excel_file)
    assert('Variant' in xlsx.get_sheet_names())
    variants = xlsx.get_sheet_by_name('Variant')
    belowHeader = False
    for row in variants.rows:
        if row[0].value == "##Local ID":
            headerRow = translateHeader(row)
            tsv.write("\t".join(headerRow) + "\n")
        elif row[0].value == "Please start your submission in the next row.":
            belowHeader = True
        elif belowHeader:
            contentsRow = translateContents(row)
            tsv.write("\t".join(contentsRow) + "\n")
        
def translateHeader(excelHeader):
    """
    Parse the header line.  Perform specific column renames.  Chop off the last column,
    return the rest.
    """
    translatedHeader = list()
    for cell in excelHeader:
        contents = cell.value
        contents = re.sub("Alternate designations", "BIC Nomenclature", contents)
        contents = re.sub("Official allele name", "Abbrev AA change", contents)
        contents = re.sub("'", "", contents)
        translatedHeader.append(contents)
    return(translatedHeader[:-1])

def translateContents(excelContents):
    """
    Parse the conents.  Chop off the last column and return the rest
    """
    translatedContents = list()
    for cell in excelContents:
        if cell.value is None:
            translatedContents.append("")
        else:
            thisValue = str(cell.value)
            translatedContents.append(thisValue)
    return(translatedContents[:-1])

if __name__ == "__main__":
    main()
