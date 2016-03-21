#!/usr/bin/env python
"""
this scripts takes the enigma variant list and merge vcf files in a folder into
the exisitng enigma variants:
"""
import argparse
import vcf
import subprocess
import tempfile
import shutil
from StringIO import StringIO
from copy import deepcopy
from pprint import pprint
import pickle


BRCA1 = {"hg38": {"start": 43000000,
                  "sequence": open("../resources/brca1_hg38.txt", "r").read()},
         "hg19": {"start": 41100000,
                  "sequence": open("../resources/brca1_hg19.txt", "r").read()}}
BRCA2 = {"hg38": {"start": 32300000,
                  "sequence": open("../resources/brca2_hg38.txt", "r").read()},
         "hg19": {"start": 32800000,
                  "sequence": open("../resources/brca2_hg19.txt", "r").read()}}
  
#GENOMIC VERSION:
VERSION = "hg38" # equivalent to GRCh38


# files needed for string comparison

#key value pair dictionaries of all extra fields in various databases to add
GENOME1K_FIELDS = {"Allele_frequency":"AF",
                   "EAS_Allele_frequency":"EAS_AF",
                   "EUR_Allele_frequency":"EUR_AF",
                   "AFR_Allele_frequency":"AFR_AF",
                   "AMR_Allele_frequency":"AMR_AF",
                   "SAS_Allele_frequency":"SAS_AF"}
CLINVAR_FIELDS = {"HGVS": "HGVS",
                  "Submitter":"Submitter",
                  "Clinical_Significance":"ClinicalSignificance",
                  "Date_Last_Updated":"DateLastUpdated",
                  "SCV":"SCV",
                  "Allele_Origin":"Origin",
                  "Protein":"Protein",
                  "Method":"Method"}
LOVD_FIELDS = {"Origin_of_variant": "genetic_origin",
               "Variant_frequency": "frequency",
               "Variant_haplotype": "haplotype",
               "Functional_analysis_result": "functionalanalysis_result",
               "Functional_analysis_technique": "functionalanalysis_technique",
               "HGVS_cDNA": "dna_change",
               "HGVS_protein": "protein_change"}
EX_LOVD_FIELDS = {"Combined_prior_probablility": "combined_prior_p",
                  "Segregation_LR": "segregation_lr",
                  "Sum_family_LR": "sum_family_lr",
                  "Co_occurrence_LR": "co_occurrence_lr",
                  "Missense_analysis_prior_probability": "missense_analysis_prior_p",
                  "Posterior_probability": "posterior_p",
                  "IARC_class":"iarc_class",
                  "BIC_Nomenclature": "bic_dna_change",
                  "Literature_source":"observational_reference",
                  "HGVS_cDNA": "dna_change",
                  "HGVS_protein": "protein_change"}
BIC_FIELDS = {"Clinical_classification": "Category",
              "Number_of_family_member_carrying_mutation": "Number_Reported",
              "Patient_nationality": "Nationality",
              "Germline_or_Somatic": "G_or_S",
              "Mutation_type": "Mutation_Type",
              "BIC_Designation": "Designation",
              "Clinical_importance": "Clinically_Importance",
              "Ethnicity": "Ethnicity",
              "Literature_citation": "Reference"}
ESP_FIELDS = {"polyPhen2_result": "PH",
              "Minor_allele_frequency":"MAF"}
EXAC_FIELDS = {"Allele_frequency": "AF"}

FIELD_DICT = {"1000_Genomes": GENOME1K_FIELDS,
               "ClinVar": CLINVAR_FIELDS,
               "LOVD": LOVD_FIELDS,
               "exLOVD": EX_LOVD_FIELDS,
               "ExAC": EXAC_FIELDS,
               "ESP": ESP_FIELDS,
               "BIC": BIC_FIELDS}

ENIGMA_FILE = "enigma_variants_GRCh38_2-27-2016.tsv"
GENOME1K_FILE = "10k_genome.brca.sorted.hg38.vcf"
CLINVAR_FILE = "ClinVarBrca.vcf"
LOVD_FILE = "sharedLOVD_brca12.sorted.hg38.vcf"
EX_LOVD_FILE = "exLOVD_brca12.sorted.hg38.vcf"
BIC_FILE = "bic_brca12.sorted.hg38.vcf"
EXAC_FILE = "exac_BRCA12.sorted.hg38.vcf"
ESP_FILE = "esp.brca.vcf"


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input VCF directory",
                    default="/hive/groups/cgl/brca/release1.0/pipeline_input/")
parser.add_argument("-o", "--output", 
                    default="/hive/groups/cgl/brca/release1.0/merged.csv")
parser.add_argument("-e", "--ev", 
                    default="/hive/groups/cgl/brca/release1.0/equivalent_variants.pickledumps")
parser.add_argument("-w", "--wrong_genome", 
                    help="Directory for data with  wrong genomic coordinates",
                    default="/hive/groups/cgl/brca/release1.0/vcf_wrong_genome_coordinate/")

ARGS = parser.parse_args()




def main():
    tmp_dir = tempfile.mkdtemp()
    try:
        source_dict, columns, variants = preprocessing(tmp_dir)
        print "------------merging different dataset------------------------------"
        for source_name, file in source_dict.iteritems():
            (columns, variants) = add_new_source(columns, variants, source_name, 
                                                 file, FIELD_DICT[source_name])
        print "------------string comparison merge-------------------------------"
        variants = string_comparison_merge(variants) 
        write_new_csv(ARGS.output, columns, variants)
        print ARGS.input
        print ARGS.output
        print ARGS.ev
        print ARGS.wrong_genome
    finally:
        shutil.rmtree(tmp_dir)

def string_comparison_merge(variants):
    # make sure the input genomic coordinate strings are already unique strings
    assert (len(variants.keys()) == len(set(variants.keys())))
    #equivalence = find_equivalent_variant(variants.keys())
    #with open(ARGS.ev, "w") as f:
    #    f.write(pickle.dumps(equivalence))
    #f.close()
    equivalence = pickle.loads(open(ARGS.ev, "r").read())
    for equivalent_v in equivalence:
        merged_row = []
        for each_v in equivalent_v:
            if len(merged_row) == 0:
                merged_row = variants[each_v]
                variants.pop(each_v)
            else:
                for index, row1 in enumerate(merged_row):
                    row2 = variants[each_v][index]
                    if row1 == "-" and row2 != "-":
                        merged_row[index] = row2
                    elif row1 != "-" and row2 == "-":
                        merged_row[index] = row1
                    else:
                        if row1 != row2 and (row2 not in row1.split("|")):
                            merged_row[index] = row1 + "|" + row2
                variants.pop(each_v)
        variants["|".join(list(equivalent_v))] = merged_row
    return variants

def find_equivalent_variant(genome_coors):
    uniq_variants = {}
    for i, v in enumerate(genome_coors):
        variant_exist = False
        for existing_v in uniq_variants:
            if v == existing_v:
                continue
            else:
                v1 = v.replace("-", "").replace("chr", "").replace(">", ":")
                v2 = existing_v.replace("-", "").replace("chr", "").replace(">", ":")
                if variant_equal(v1.split(":"), v2.split(":")):
                    variant_exist = True
                    uniq_variants[existing_v].add(v)
        if not variant_exist:
            uniq_variants[v] = set([v])
    equivalent_variants = [] 
    for value in uniq_variants.values():
        if len(value) > 1:
            equivalent_variants.append(value)
    return equivalent_variants

def preprocessing(tmp_dir):
    # Preprocessing variants:
    source_dict = {"1000_Genomes": GENOME1K_FILE + "for_pipeline",
                   "ClinVar": CLINVAR_FILE,
                   "LOVD": LOVD_FILE,
                   "exLOVD": EX_LOVD_FILE,
                   "ExAC": EXAC_FILE,
                   "ESP": ESP_FILE,
                   "BIC": BIC_FILE,
                   }    
    print "\n" + ARGS.input + ":"
    print "ENIGMA: {0}".format(ENIGMA_FILE)
    for source_name, file_name in source_dict.iteritems():
        print source_name, ":", file_name
    print "------------preprocessing--------------------------------"
    print "remove sample columns and two erroneous rows from 1000 Genome file"
    f_1000G = open(ARGS.input + GENOME1K_FILE + "for_pipeline", "w")
    subprocess.call(
       ["bash", "1000g_preprocess.sh", ARGS.input + GENOME1K_FILE], stdout=f_1000G)
    
    print "-------check if genomic coordinates are correct----------"
    (columns, variants) = save_enigma_to_dict(ARGS.input + ENIGMA_FILE)
    for source_name, file_name in source_dict.iteritems():
        f = open(ARGS.input + file_name, "r")
        f_wrong = open(ARGS.wrong_genome + source_name + "_wrong_genome_coor.vcf", "w")
        f_right = open(tmp_dir + "/right" + source_name, "w")
        vcf_reader = vcf.Reader(f, strict_whitespace=True)
        vcf_wrong_writer = vcf.Writer(f_wrong, vcf_reader)
        vcf_right_writer = vcf.Writer(f_right, vcf_reader)
        n_wrong, n_total = 0, 0
        for record in vcf_reader:
            v = [record.CHROM, record.POS, record.REF, record.ALT]
            if not ref_correct(v):
                vcf_wrong_writer.write_record(record)
                n_wrong += 1
            else:
                vcf_right_writer.write_record(record)
            n_total += 1
        f_right.close()
        f_wrong.close()
        print "in {0}, wrong: {1}, total: {2}".format(source_name, n_wrong, n_total) 
    print "variants with wrong genomic coordates are saved to:", ARGS.wrong_genome
    print "---------------------------------------------------------"
    
    for source_name, file_name in source_dict.iteritems():
        print "convert to one variant per line in ", source_name
        f_in = open(tmp_dir + "/right" + source_name, "r")
        f_out = open(tmp_dir + "/" + source_name + ".vcf", "w")
        one_variant_transform(f_in, f_out)
        print "merge repetitive variants within ", source_name
        f_in = open(tmp_dir + "/" + source_name + ".vcf", "r")
        f_out = open(tmp_dir + "/" + source_name + "ready.vcf", "w")
        repeat_merging(f_in, f_out)
        source_dict[source_name] = f_out.name
    return source_dict, columns, variants

def repeat_merging(f_in, f_out):
    """takes a vcf file, collapses repetitive variant rows and write out
        to a new vcf file (without header)"""
    vcf_reader = vcf.Reader(f_in, strict_whitespace=True)
    variant_dict = {}
    num_repeats = 0
    for record in vcf_reader:
        genome_coor = "chr{0}:{1}:{2}>{3}".format(
            record.CHROM, str(record.POS), record.REF, record.ALT[0])
        if genome_coor not in variant_dict.keys():
            variant_dict[genome_coor] = deepcopy(record)
        else:
            num_repeats += 1
            for key in record.INFO:        
                if key not in variant_dict[genome_coor].INFO.keys():
                    variant_dict[genome_coor].INFO[key] = deepcopy(record.INFO[key])
                else:
                    new_value = deepcopy(record.INFO[key])
                    old_value = deepcopy(variant_dict[genome_coor].INFO[key])

                    if type(new_value) != list:
                        new_value = [new_value]
                    if type(old_value) != list:
                        old_value = [old_value]
                    if new_value  == old_value:
                        continue
                    else:
                        merged_value = list(set(new_value + old_value))
                        variant_dict[genome_coor].INFO[key] = deepcopy(merged_value)
    print "number of repeat records: ", num_repeats, "\n"
    write_to_vcf(f_out, variant_dict)

def write_to_vcf(f_out, v_dict):
    for record in v_dict.values():
        if record.QUAL == None:
            record.QUAL = "."
        if record.FILTER == None:
            record.FILTER = "."

        items = [record.CHROM, str(record.POS), str(record.ID), record.REF, 
                str(record.ALT[0]), record.QUAL, record.FILTER]
        infos = []
        for key in record.INFO:
            this_info = record.INFO[key] 
            if type(this_info) == list:
                this_info = [str(x) for x in this_info]
                infos.append(key + "=" + "|".join(this_info))
            else:
                infos.append(key + "=" + str(this_info))
        items.append(";".join(infos))
        new_line = "\t".join([str(i) for i in items])
        f_out.write(new_line + "\n")

def one_variant_transform(f_in, f_out):
    """takes a vcf file, read each row, if the ALT field contains more than 
       one item, create multiple variant row based on that row, writes new vcf"""
    vcf_reader = vcf.Reader(f_in, strict_whitespace=True)
    vcf_writer = vcf.Writer(f_out, vcf_reader)
    for record in vcf_reader:
        n = len(record.ALT)
        if n == 1:
            vcf_writer.write_record(record)
        else:
            for i in range(n):
                new_record = deepcopy(record)
                new_record.ALT = [deepcopy(record.ALT[i])]
                for key in record.INFO.keys():
                    value = deepcopy(record.INFO[key])
                    if type(value) == list and len(value) == n:
                        new_record.INFO[key] = [value[i]]
                vcf_writer.write_record(new_record)

def write_new_csv(filename, columns, variants):
    merged_file = open(filename, "w")
    merged_file.write(",".join(columns)+"\n")
    for variant in variants.values():
        if len(variant) != len(columns):
            raise Exception("mismatching number of columns in head and row")
        merged_file.write(",".join(variant)+"\n")

def add_new_source(columns, variants, source, source_file, source_dict):
    print "adding {0} into merged file.....".format(source)
    old_column_num = len(columns)
    for column_title in source_dict.keys():
        columns.append(column_title+"({0})".format(source))

    vcf_reader = vcf.Reader(open(source_file, 'r'), strict_whitespace=True)
    overlap = 0
    variants_num = 0
    for record in vcf_reader:
        variants_num += 1
        genome_coor = ("chr" + str(record.CHROM) + ":" + str(record.POS) + ":" +
                       record.REF + ">" + str(record.ALT[0]))
        if genome_coor in variants.keys():
            overlap += 1
            variants[genome_coor][0] += "|{0}".format(source)
        else: 
            variants[genome_coor] = ['-'] * old_column_num
            variants[genome_coor][0] = source
            variants[genome_coor][2] = genome_coor

        for value in source_dict.values():
            try:
                variants[genome_coor].append(str(record.INFO[value]))
            except KeyError:
                if source == "BIC":
                    variants[genome_coor].append("-")
                else:
                    raise Exception("uncaught weirdness")

    # for those enigma record that doesn't have a hit with new genome coordinate
    # add extra cells of "-" to the end of old record
    for value in variants.values():
        if len(value) != len(columns):
            value += ["-"] * len(source_dict)

    print "number of variants in " + source + " is ", variants_num
    print "overlap with previous dataset: ", overlap
    print "number of total variants with the addition of " + source + " is: ", len(variants), "\n"

    for value in variants.values():
        if len(value) != len(columns):
            raise Exception("mismatching number of columns in head and row")
    return (columns, variants)


def save_enigma_to_dict(path):
    enigma_file = open(path, "r")
    variants = dict()
    columns = ""
    line_num = 0
    f_wrong = open(ARGS.wrong_genome + "ENIGMA_wrong_genome.txt", "w")
    n_wrong, n_total = 0, 0
    for line in enigma_file:
        line_num += 1
        line = line.replace(",", "|")
        if line_num == 1:
            columns = line.strip().split("\t")
            columns = [c + "(ENIGMA)" for c in columns if c != "Genomic_Coordinate"]
            columns.insert(0, "Source")
            columns.insert(2, "Genomic_Coordinate")
            f_wrong.write(line)
        else:
            items = line.strip().split("\t")
            items.insert(0, "ENIGMA")
            v = items[2].replace("-", "").replace("chr", "").replace(">", ":")
            if ref_correct(v.split(":")):
                variants[items[2]] = items
            else:
                n_wrong += 1
                f_wrong.write(line)
            n_total += 1
    f_wrong.close()
    print "in ENIGMA, wrong: {0}, total: {1}".format(n_wrong, n_total)
    return (columns, variants)


def variant_equal(v1, v2, version="hg38"):
    " return (edited1, edited2) "
    if v1 == v2:
        return True

    chr1, pos1, ref1, alt1 = v1
    chr2, pos2, ref2, alt2 = v2
    pos1 = int(pos1)
    pos2 = int(pos2)

    if chr1 != chr2:
        return False
    if (len(ref1) - len(alt1)) != (len(ref2) - len(alt2)):
        return False

    # if len(ref2)>100 or len(ref1)>100:
    #     return False

    # make sure that v1 is upstream of v2
    if pos1 > pos2:
        return variant_equal(v2, v1)

    # lift coordinates and make everything 0-based
    if chr1 == "13":
        seq = BRCA2[version]["sequence"]
        pos1 = pos1 - 1 - BRCA2[version]["start"]
        pos2 = pos2 - 1 - BRCA2[version]["start"]
    elif chr1 == "17":
        seq = BRCA1[version]["sequence"]
        pos1 = pos1 - 1 - BRCA1[version]["start"]
        pos2 = pos2 - 1 - BRCA1[version]["start"]
    else:
        assert(False)

    # replace vcf ref string with alt string
    edited_v1 = seq[0:pos1]+alt1+seq[pos1+len(ref1):]
    edited_v2 = seq[0:pos2]+alt2+seq[pos2+len(ref2):]

    return edited_v1 == edited_v2

def ref_correct(v, version="hg38"): 
    chr, pos, ref, alt = v 
    pos = int(pos) 
    if chr == "13": 
        seq = BRCA2[version]["sequence"] 
        pos = pos - 1 - BRCA2[version]["start"]
    elif chr == "17":
        seq = BRCA1[version]["sequence"]
        pos = pos - 1 - BRCA1[version]["start"]
    else:
        assert(False)

    genomeRef = seq[pos:pos+len(ref)].upper()
    if len(ref) != 0 and len(genomeRef)==0:
        print v
        raise Exception("ref not inside BRCA1 or BRCA2")
    if (genomeRef != ref):
        return False
    else:
        return True


if __name__ == "__main__":
    main()

