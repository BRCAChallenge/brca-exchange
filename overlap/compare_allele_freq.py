import overlap_string_comparison as overlap2
import venn_diagram as v

"""
this code compares the allele frequency of the 909 common variants between ClinVar, UMD and LOVD
structure of variant_dic:

 {"variant1":
    {"ClinVar":[ClinVar_info1, ClinVar_info2...]
     "UMD": [UMD_info1, UMD_info2, ...]
     "LOVD": [LOVD_info1, LOVD_info2, ...]
     }
  "variant2":
    {"ClinVar":[ ]
     "UMD": [ ]
     "LOVD" : [ ]
    }
 }

"""

def main():
    varF = open("909_common_variants", "r")
    clinvarV = overlap2.get_unique_variants("allele_freq_vcfs/ClinVar.allele_freq.vcf")
    lovdV = overlap2.get_unique_variants("allele_freq_vcfs/LOVD.allele_freq.vcf")
    umdV = overlap2.get_unique_variants("allele_freq_vcfs/UMD.allele_freq.vcf")

    print len(clinvarV)
    print len(lovdV)
    print len(umdV)
    n1 = 0
    n2 = 0
    n3 = 0
    sets = {}

    sets["A"] = []
    sets["B"] = []
    sets["C"] = []


    for line in varF:
         this_variant = line.strip().split(".")
         if overlap2.check_variant_exist(this_variant, clinvarV):
             sets["A"].append(this_variant)
         if overlap2.check_variant_exist(this_variant, lovdV):
             sets["B"].append(this_variant)
         if overlap2.check_variant_exist(this_variant, umdV):
             sets["C"].append(this_variant)
    sets["AB"] = v.get_overlap(sets["A"], sets["B"])
    sets["AC"] = v.get_overlap(sets["A"], sets["C"])
    sets["BC"] = v.get_overlap(sets["B"], sets["C"])
    sets["ABC"] = v.get_overlap(sets["AB"], sets["C"])

    for key in sets:
        print key, len(sets[key])






if __name__ == "__main__":
    main()
    