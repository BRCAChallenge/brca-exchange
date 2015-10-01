"""
this scripts takes the enigma variant list and merge vcf files in a folder into
the exisitng enigma variants:
"""
import vcf

"key value pair dictionaries of all extra fields in various databases to add"
GENOME1K = {"Allele_frequency":"AF", "EAS_Allele_frequency":"EAS_AF", "EUR_Allele_frequency":"EUR_AF",
              "AFR_Allele_frequency":"AFR_AF", "AMR_Allele_frequency":"AMR_AF", "SAS_Allele_frequency":"SAS_AF"}
CLINVAR = {"Allele_origin":"CLNORIGIN", 
           #"Nonsense_mutation": "NSN", "Missense_mutation":"NSM", "Frameshit_mutation":"NSF",
           "Variant_clinical_significance":"CLNSIG"} 




def main():
    (columns, variants) = save_enigma_to_dict("data/enigma_variants_9-29-2015.tsv")
    print "number of variants in enigma", len(variants)
#    (columns, variants) = add_new_source(columns, variants, "1000_Genomes", 
#                                         'data/allVcf/1000_genomes.brca.vcf', GENOME1K) 
    print "number of variants in enigma + 1000 genomes", len(variants)
    
    (columns, variants) = add_new_source(columns, variants, "ClinVar", 
                                         'data/allVcf/brca.clinvar.vcf', CLINVAR) 
    print "number of variants in enigma + Clinvar", len(variants)
    write_new_tsv("data/merge/enigma-1000genomes-clinvar.tsv", columns, variants)       
 
def write_new_tsv(filename, columns, variants):
    merged_file = open(filename, "w")
    merged_file.write("\t".join(columns)+"\n")
    for variant in variants.values():
        merged_file.write("\t".join(variant)+"\n")

def add_new_source(columns, variants, source, source_file, source_dict):
    num_columns = len(source_dict)
    for column_title in source_dict.keys():
        columns.append(column_title)
    for key in variants.keys():
        variants[key] += ["-"] * num_columns

    vcf_reader = vcf.Reader(open(source_file, 'r'))
    overlap = 0
    variants_num = 0
    for record in vcf_reader:
        for i in range(len(record.ALT)):
            variants_num += 1
            genome_coor = (record.CHROM + ":" + str(record.POS) + ":" + 
                           record.REF + ">" + str(record.ALT[i]))
            if genome_coor in variants.keys():
                overlap += 1
                ### TODO: add actual code to make this part work
                variants[genome_coor][0] += ",ClinVar"
            else: 
                variants[genome_coor] = ['-'] * len(columns)
                variants[genome_coor][0] = source 
                variants[genome_coor][2] = genome_coor
            
            for value in source_dict.values():
                if len(record.INFO[value]) == len(record.ALT):
                    variants[genome_coor].append(str(record.INFO[value][i]))
                else:
                    variants[genome_coor].append(str(record.INFO[value][0]))

    print "number of variants in " + source + " is ", variants_num
    print "overlap: ", overlap
    return (columns, variants) 


def save_enigma_to_dict(path):
    enigma_file = open(path, "r")
    variants = dict()
    columns = ""
    line_num = 0
    for line in enigma_file:
        line_num += 1
        if line_num == 1:
            columns = line.strip().split("\t")
            columns.insert(0,"Source")
        else:
            items = line.strip().split("\t")
            items.insert(0,"ENIGMA")
            variants[items[2]] = items
    return (columns, variants)

if __name__=="__main__":
    main()

