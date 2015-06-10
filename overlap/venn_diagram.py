from matplotlib_venn import venn3
from matplotlib_venn import venn2
from matplotlib import pyplot as plt
import sys

import overlap_string_comparison as overlap2
import overlap_direct_comparison as overlap1



"""
    venn3() takes a list of 7 numbers:
    venn[0] -> number exclusively in A
    venn[1] -> number exclusively in B
    venn[2] -> overlap between A and B but not C
    venn[3] -> number exclusively in C
    venn[4] -> overlap between A and C but not B
    venn[5] -> overlap between B and C but not A
    venn[6] -> overlap between A, B and C
    plot:  venn3(subsets=venn, set_labels = ("A", "B", "C")
"""

PATH = "/hive/groups/cgl/brca/phase1/data/cutoff_vcf/"


def main():
    A = "ClinVar"
    B = "LOVD"
    C = "UMD" # optional, leave empty "" if only wanting two circles

    # sets = get_string_comparison_set(A, B, C, PATH)
    sets = get_direct_comparison_set(A, B, C, PATH)


    if C == "": # draw two-circle venn diagram
        draw_venn2(A, B, sets)
    else:
        draw_venn3(A, B, C, sets)


def get_string_comparison_set(A, B, C, path):
    sets = {}
    sets["A"] = overlap2.get_unique_variants(path + A + ".brca.no_header.cutoff.vcf")
    sets["B"] = overlap2.get_unique_variants(path + B + ".brca.no_header.cutoff.vcf")
    sets["AB"] = get_overlap(sets["A"], sets["B"])
    if C != "":
        sets["C"] = overlap2.get_unique_variants(path + C + ".brca.no_header.cutoff.vcf")
        sets["BC"] = get_overlap(sets["B"], sets["C"])
        sets["AC"] = get_overlap(sets["A"], sets["C"])
        sets["ABC"] = get_overlap(sets["AB"], sets["C"])
    return sets

def get_direct_comparison_set(A, B, C, path):
    sets = {}
    sets["A"] = set(overlap1.get_variant_list(path + A + ".brca.no_header.cutoff.vcf"))
    sets["B"] = set(overlap1.get_variant_list(path + B + ".brca.no_header.cutoff.vcf"))
    sets["AB"] = sets["A"].intersection(sets["B"])
    if C != "":
        sets["C"] = set(overlap1.get_variant_list(path + C + ".brca.no_header.cutoff.vcf"))
        sets["AC"] = sets["A"].intersection(sets["C"])
        sets["BC"] = sets["B"].intersection(sets["C"])
        sets["ABC"] = sets["AB"].intersection(sets["C"])
    return sets


def draw_venn3(A, B, C, sets):
    venn = [0]*7
    venn[2] = len(sets["AB"]) - len(sets["ABC"])
    venn[4] = len(sets["AC"]) - len(sets["ABC"])
    venn[5] = len(sets["BC"]) - len(sets["ABC"])
    venn[6] = len(sets["ABC"])
    venn[0] = len(sets["A"]) - venn[2] - venn[4] - venn[6]
    venn[1] = len(sets["B"]) - venn[2] - venn[5] - venn[6]
    venn[3] = len(sets["C"]) - venn[4] - venn[5] - venn[6]

    labelA = A + " (" + str(len(sets["A"])) + ")"
    labelB = B + " (" + str(len(sets["B"])) + ")"
    labelC = C + " (" + str(len(sets["C"])) + ")"
    venn3(subsets=venn, set_labels = (labelA, labelB, labelC))
    plt.show()

def draw_venn2(A, B, sets):
    venn = [0]*3
    venn[2] = len(sets["AB"])
    venn[0] = len(sets["A"]) - len(sets["AB"])
    venn[1] = len(sets["B"]) - len(sets["AB"])
    labelA = A + " (" + str(len(sets["A"])) + ")"
    labelB = B + " (" + str(len(sets["B"])) + ")"
    print venn
    venn2(subsets=venn, set_labels=(labelA, labelB))
    plt.show()


def get_overlap(db1, db2):
    variants = []
    for variant in db1:
        if overlap2.check_variant_exist(variant, db2):
            variants.append(variant)
    return variants


if __name__ == "__main__" :
    sys.exit(main())

