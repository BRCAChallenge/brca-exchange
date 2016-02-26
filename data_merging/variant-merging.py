"""
this scripts takes the enigma variant list and merge vcf files in a folder into
the exisitng enigma variants:
"""
import vcf

#key value pair dictionaries of all extra fields in various databases to add
GENOME1K = {"Allele_frequency(1000_Genomes)":"AF",
            "EAS_Allele_frequency(1000_Genomes)":"EAS_AF",
            "EUR_Allele_frequency(1000_Genomes)":"EUR_AF",
            "AFR_Allele_frequency(1000_Genomes)":"AFR_AF",
            "AMR_Allele_frequency(1000_Genomes)":"AMR_AF",
            "SAS_Allele_frequency(1000_Genomes)":"SAS_AF"}
CLINVAR = {"Allele_origin(ClinVar)":"CLNORIGIN",
           "Variant_clinical_significance(ClinVar)":"CLNSIG"}
LOVD = {"Origin_of_variant(LOVD)": "genetic_origin",
        "Variant_frequency(LOVD)": "frequency",
        "Variant_haplotype(LOVD)": "haplotype",
        "Variant_affecting_protein(LOVD)": "effect",
        "HGVS_cDNA(LOVD)": "dna_change",
        "HGVS_genomic(LOVD)": "dna_change_genomic",
        "HGVS_protein(LOVD)": "protein_change"}
EXAC = {"Allele_frequency(ExAC)": "AF",
        "VEP_Gene(ExAC)": "CSQ_Gene",
        "VEP_Consequence(ExAC)":"CSQ_Consequence",
        "VEP_HGVSc(ExAC)":"CSQ_HGVSc",
        "VEP_HGVSp(ExAC)":"CSQ_HGVSp"}
EX_LOVD = {"Exon_number(exLOVD)":"exon",
           "HGVS_cDNA(exLOVD)":"dna_change",
           "BIC(exLOVD)":"dna_change_bic",
           "HGVS_protein(exLOVD)":"protein_change",
           "IARC_class(exLOVD)":"iarc_class",
           "Literature_source(exLOVD)":"observational_reference"}
BIC = {"Clinical_classification(BIC)":"Category",
       "Number_of_family_member_carrying_mutation(BIC)":"Number_Reported",
       "Exon_number(BIC)":"Exon",
       "Patient_nationality(BIC)":"Nationality",
       "Germline_or_Somatic(BIC)":"G_or_S",
       "Mutation_type(BIC)":"Mutation_Type",
       "BIC_Designation(BIC)":"Designation",
       "Clinical_importance(BIC)":"Clinically_Importance",
       "Ethnicity(BIC)":"Ethnicity",
       "HGVS_cDNA(BIC)":"HGVS_cDNA",
       "Literature_citation(BIC)":"Reference",
       "HGVS_genomic(BIC)":"HGVS_Genomic_(hg19)",
       "HGVS_protein(BIC)":"HGVS_Protein"}

ENIGMA_FILE = "../data/enigma_variants_9-29-2015.tsv"
GENOME1K_FILE = "../data/allVcf/no_repeats/1000_genomes.brca.no_sample.ovpr.no_repeats.vcf"
CLINVAR_FILE = "../data/allVcf/no_repeats/clinvar.brca.ovpr.no_repeats.vcf"
LOVD_FILE = "../data/allVcf/no_repeats/lovd.brca.ovpr.no_repeats.vcf"
EX_LOVD_FILE = "../data/allVcf/no_repeats/ex_lovd.brca.ovpr.no_repeats.vcf"
BIC_FILE = "../data/allVcf/no_repeats/bic.brca.no_repeats.vcf"
EXAC_FILE = "../data/allVcf/no_repeats/exac.brca.ovpr.no_repeats.vcf"

SOURCE_DICT = {"1000_Genomes": [GENOME1K_FILE, GENOME1K],
               "ClinVar": [CLINVAR_FILE, CLINVAR],
               "LOVD": [LOVD_FILE, LOVD],
               "ExAC": [EXAC_FILE, EXAC],
               "BIC": [BIC_FILE, BIC]}


def main():
    (columns, variants) = save_enigma_to_dict(ENIGMA_FILE)

    for source, value in SOURCE_DICT.iteritems():
        (columns, variants) = add_new_source(columns, variants, source,
                                             value[0], value[1])

#    (columns, variants) = add_new_source(columns, variants, "1000_Genomes",
#                                          GENOME1K_FILE, GENOME1K)
#
#    (columns, variants) = add_new_source(columns, variants, "ClinVar",
#            CLINVAR_FILE, CLINVAR)
#
#    (columns, variants) = add_new_source(columns, variants, "LOVD",
#            LOVD_FILE, LOVD)
#
#    (columns, variants) = add_new_source(columns, variants, "ExAC",
#            EXAC_FILE, EXAC)
#
#    (columns, variants) = add_new_source(columns, variants, "exLOVD",
#            EX_LOVD_FILE, EX_LOVD)
#
#    (columns, variants) = add_new_source(columns, variants, "BIC",
#            BIC_FILE, BIC)

    write_new_tsv("../data/merge/merged_new.tsv", columns, variants)

def write_new_tsv(filename, columns, variants):
    merged_file = open(filename, "w")
    merged_file.write("\t".join(columns)+"\n")
    for variant in variants.values():
        if len(variant) != len(columns):
            raise Exception("mismatching number of columns in head and row")
        merged_file.write("\t".join(variant)+"\n")

def add_new_source(columns, variants, source, source_file, source_dict):
    old_column_num = len(columns)
    for column_title in source_dict.keys():
        columns.append(column_title)

    vcf_reader = vcf.Reader(open(source_file, 'r'), strict_whitespace=True)
    overlap = 0
    variants_num = 0
    for record in vcf_reader:
        variants_num += 1
        genome_coor = ("chr" + str(record.CHROM) + ":" + str(record.POS) + ":" +
                       record.REF + ">" + str(record.ALT[0]))
        if genome_coor in variants.keys():
            overlap += 1
            variants[genome_coor][0] += ",{0}".format(source)

        else: 
            variants[genome_coor] = ['-'] * old_column_num
            variants[genome_coor][0] = source
            variants[genome_coor][2] = genome_coor

        for value in source_dict.values():
            variants[genome_coor].append(str(record.INFO[value]))

    # for those enigma record that doesn't have a hit with new genome coordinate
    # add extra cells of "-" to the end of old record
    for value in variants.values():
        if len(value) != len(columns):
            value += ["-"] * len(source_dict)

    print "number of variants in " + source + " is ", variants_num
    print "overlap with previous dataset: ", overlap
    print "number of variants with the addition of " + source + "is: ", len(variants), "\n"


    for value in variants.values():
        if len(value) != len(columns):
            raise Exception("mismatching number of columns in head and row")

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
            columns.insert(0, "Source")
        else:
            items = line.strip().split("\t")
            items.insert(0, "ENIGMA")
            variants[items[2]] = items
    print "number of variants in enigma: ", len(variants), "\n"
    return (columns, variants)

if __name__ == "__main__":
    main()

