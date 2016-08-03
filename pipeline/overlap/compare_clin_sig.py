import overlap_string_comparison as overlap2
import venn_diagram as v

"""
this code compares the clinical significance of the 1149 common variants between ClinVar, UMD and BIC
structure of variant_dic:

 {"variant1":
    {"ClinVar":[ClinVar_info1, ClinVar_info2...]
     "BIC": [BIC_info1, BIC_info2, ...]
     "UMD": [UMD_info1, UMD_info2, ...]
     }
  "variant2":
    {"ClinVar":[ ]
     "BIC": [ ]
     "UMD" : [ ]
    }
 }

"""

#PATH = "/hive/groups/cgl/brca/phase1/data/allele_freq_vcfs/"
PATH = "/Users/charlesmarkello/clin_sig_vcfs/"

def main():
    #varF = open("909_common_variants", "r")
    varF = open("ClinVar_UMD_BIC_common_variants", "r")
    clinvarV = overlap2.get_unique_variants(PATH + "ClinVar.clin_sig.vcf")
    bicV = overlap2.get_unique_variants(PATH + "BIC.clin_sig.vcf")
    umdV = overlap2.get_unique_variants(PATH + "UMD.clin_sig.vcf")

    print len(clinvarV)
    print len(bicV)
    print len(umdV)
    sets = {}

    sets["A"] = []
    sets["B"] = []
    sets["C"] = []


    for line in varF:
         this_variant = line.strip().split(".")
         if overlap2.check_variant_exist(this_variant, clinvarV):
             sets["A"].append(this_variant)
         if overlap2.check_variant_exist(this_variant, bicV):
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
    
