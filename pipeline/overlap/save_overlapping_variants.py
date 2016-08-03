import venn_diagram as v

"""
this code writes out the result of the overlaps between database ClinVar, UMD, and LOVD to local disk

"""

PATH = "/hive/groups/cgl/brca/phase1/data/cutoff_vcf/"


def main():
    A = "ClinVar"
    B = "LOVD"
    C = "UMD"
    sets = v.get_string_comparison_set(A, B, C, PATH)
    f = open("909_common_variants", "w")
    for variant in sets["ABC"]:
        f.write("{0}.{1}.{2}.{3}\n".format(variant[0], variant[1], variant[2], variant[3]))


if __name__ == "__main__":
    main()
    