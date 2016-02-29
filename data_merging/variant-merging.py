"""
this scripts takes the enigma variant list and merge vcf files in a folder into
the exisitng enigma variants:
"""
import vcf
import subprocess
import tempfile
import shutil
from StringIO import StringIO
from copy import deepcopy


#key value pair dictionaries of all extra fields in various databases to add
GENOME1K_FIELDS = {"Allele_frequency(1000_Genomes)":"AF",
                   "EAS_Allele_frequency(1000_Genomes)":"EAS_AF",
                   "EUR_Allele_frequency(1000_Genomes)":"EUR_AF",
                   "AFR_Allele_frequency(1000_Genomes)":"AFR_AF",
                   "AMR_Allele_frequency(1000_Genomes)":"AMR_AF",
                   "SAS_Allele_frequency(1000_Genomes)":"SAS_AF"}
CLINVAR_FIELDS = {"Allele_origin(ClinVar)":"CLNORIGIN",
                  "Variant_clinical_significance(ClinVar)":"CLNSIG"}
LOVD_FIELDS = {"Origin_of_variant(LOVD)": "genetic_origin",
               "Variant_frequency(LOVD)": "frequency",
               "Variant_haplotype(LOVD)": "haplotype",
               "Variant_affecting_protein(LOVD)": "effect",
               "HGVS_cDNA(LOVD)": "dna_change",
               "HGVS_genomic(LOVD)": "dna_change_genomic",
               "HGVS_protein(LOVD)": "protein_change"}
EXAC_FIELDS = {"Allele_frequency(ExAC)": "AF"}
EX_LOVD_FIELDS = {"Exon_number(exLOVD)":"exon",
                  "HGVS_cDNA(exLOVD)":"dna_change",
                  "HGVS_protein(exLOVD)":"protein_change",
                  "IARC_class(exLOVD)":"iarc_class",
                  "Literature_source(exLOVD)":"observational_reference"}
FIELD_DICT = {"1000_Genomes": GENOME1K_FIELDS,
               #"ClinVar": CLINVAR_FIELDS,
               "LOVD": LOVD_FIELDS,
               "exLOVD": EX_LOVD_FIELDS,
               "ExAC": EXAC_FIELDS,
               #"ESP": ESP_FILE,
              }

PIPELINE_INPUT = "/hive/groups/cgl/brca/release1.0/pipeline_input/"
ENIGMA_FILE = PIPELINE_INPUT + "enigma_variants_GRCh38_2-27-2016.tsv"
GENOME1K_FILE = PIPELINE_INPUT + "1000G_brca.sorted.hg38.vcf"
#CLINVAR_FILE = PIPELINE_INPUT + "ClinVarBrca.vcf"
LOVD_FILE = PIPELINE_INPUT + "sharedLOVD_brca12.sorted.hg38.vcf"
EX_LOVD_FILE = PIPELINE_INPUT + "exLOVD_brca12.sorted.hg38.vcf"
EXAC_FILE = PIPELINE_INPUT + "exac_BRCA12.sorted.hg38.vcf"




def main():
    tmp_dir = tempfile.mkdtemp()
    try:
        source_dict = preprocessing(tmp_dir)
        print "------------merging------------------------------------"
        (columns, variants) = save_enigma_to_dict(ENIGMA_FILE)
        for source_name, file in source_dict.iteritems():
            (columns, variants) = add_new_source(columns, variants, source_name, 
                                                 file, FIELD_DICT[source_name])
        write_new_tsv(PIPELINE_INPUT + "merged.tsv", columns, variants)
    finally:
        shutil.rmtree(tmp_dir)


def preprocessing(tmp_dir):
    # Preprocessing variants:
    print "------------preprocessing--------------------------------"
    tmp_1000g = tmp_dir + "/1000G.vcf"
    f_1000G = open(tmp_1000g, "w")
    print "remove sample columns and two erroneous rows from 1000 Genome file"
    (subprocess.call(["bash", "preprocess.sh", GENOME1K_FILE], stdout=f_1000G))
    source_dict = {"1000_Genomes": tmp_1000g,
                   #"ClinVar": CLINVAR_FILE,
                   "LOVD": LOVD_FILE,
                   "exLOVD": EX_LOVD_FILE,
                   "ExAC": EXAC_FILE,
                   #"ESP": ESP_FILE,
                   }    
    for source_name, file_name in source_dict.iteritems():
        print "convert to one variant per line in ", source_name
        f_in = open(file_name, "r")
        f_out = open(tmp_dir + "/" + source_name + ".vcf", "w")
        one_variant_transform(f_in, f_out)
        print "merge repetitive variants within ", source_name
        f_in = open(tmp_dir + "/" + source_name + ".vcf", "r")
        f_out = open(tmp_dir + "/" + source_name + "ready.vcf", "w")
        repeat_merging(f_in, f_out)
        source_dict[source_name] = f_out.name
    return source_dict

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

def write_new_tsv(filename, columns, variants):
    merged_file = open(filename, "w")
    merged_file.write("\t".join(columns)+"\n")
    for variant in variants.values():
        if len(variant) != len(columns):
            raise Exception("mismatching number of columns in head and row")
        merged_file.write("\t".join(variant)+"\n")

def add_new_source(columns, variants, source, source_file, source_dict):
    print "adding {0} into merged file.....".format(source)
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

