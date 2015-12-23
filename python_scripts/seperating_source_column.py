"""
seperating the Variant_Source column merged_v4.tsv into six columns:
Variant_in_ENIGMA
Variant_in_ClinVar
Variant_in_1000_Genomes
Variant_in_ExAC
Variant_in_LOVD
Variant_in_BIC
"""

COLUMNS = ["Variant_in_ENIGMA",
           "Variant_in_ClinVar",
           "Variant_in_1000_Genomes",
           "Variant_in_ExAC",
           "Variant_in_LOVD",
           "Variant_in_BIC"]

SOURCES = ["ENIGMA", "ClinVar", "1000_Genomes", "ExAC", "LOVD", "BIC"]


def main():
    f_in = open('/Users/Molly/Desktop/BRCA Research/data/merged_v4.tsv', "r")
    f_out = open('/Users/Molly/Desktop/BRCA Research/data/merged_v5.tsv', "w")
    line_num = 1
    for line in f_in:
        items = line.strip().split("\t")
        if line_num == 1:
            items.pop(0)
            items = COLUMNS + items
        else:
            source = items[0]
            items.pop(0)
            items = separate_source(source) + items
        line_num += 1
        new_line = "\t".join(items) + "\n"
        f_out.write(new_line)

def separate_source(source):
    boolean_list = []
    exisiting_sources = source.split("|")
    for each_source in SOURCES:
        if each_source in exisiting_sources:
            boolean_list.append("True")
        else:
            boolean_list.append("False")
    return boolean_list


if __name__=="__main__":
    main()