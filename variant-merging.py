"""
this scripts takes the enigma variant list and merge vcf files in a folder into
the exisitng enigma variants:
"""
import vcf


def main():
    (columns, variants) = save_enigma_to_dict("data/enigma_variants_9-21-2015.tsv")
    print "number of variants in enigma", len(variants)
    (columns, variants) = add_1000_genomes(columns, variants) 
    print "number of variants in enigma + 1000 genomes", len(variants)
    print columns
    write_new_tsv("data/merge/enigma-1000genomes.tsv", columns, variants)       
 
def write_new_tsv(filename, columns, variants):
    merged_file = open(filename, "w")
    merged_file.write("\t".join(columns)+"\n")
    for variant in variants.values():
        merged_file.write("\t".join(variant)+"\n")

def add_1000_genomes(columns, variants):
    (columns, variants) = add_extra_columns(columns, variants)
    vcf_reader = vcf.Reader(open('data/allVcf/1000_genomes.brca.vcf', 'r'))
    num_1000G = 0
    overlap = 0
    for record in vcf_reader:
        for i in range(len(record.ALT)):
            num_1000G += 1
            genome_coor = (record.CHROM + ":" + str(record.POS) + ":" + 
                           record.REF + ">" + str(record.ALT[i]))
            if genome_coor in variants.keys():
                overlap += 1
                print variants[genome_coor]
            else: 
                variants[genome_coor] = ['-'] * len(columns)
                variants[genome_coor][-7] = "1000_Genomes"
                variants[genome_coor][-8] = genome_coor
            variants[genome_coor][-6] = str(record.INFO['AF'][i])
            variants[genome_coor][-5] = str(record.INFO['EAS_AF'][i])
            variants[genome_coor][-4] = str(record.INFO['EUR_AF'][i])
            variants[genome_coor][-3] = str(record.INFO['AFR_AF'][i])
            variants[genome_coor][-2] = str(record.INFO['AMR_AF'][i])
            variants[genome_coor][-1] = str(record.INFO['SAS_AF'][i])
    print "overlap between enigma and 1000 genomes: " + str(overlap)
    print "number of variants in 1000 genomes", num_1000G
    return (columns, variants) 

def add_extra_columns(columns, variants):
    columns.append("Allele_frequency")
    columns.append("EAS_Allele_frequency")
    columns.append("EUR_Allele_frequency")
    columns.append("AFR_Allele_frequency")
    columns.append("AMR_Allele_frequency")
    columns.append("SAS_Allele_frequency")
    for key in variants.keys():
        variants[key] += ["-"] * 6
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
            columns.append("Source")
        else:
            items = line.strip().split("\t")
            items.append("ENIGMA")
            variants[items[-2]] = items
    return (columns, variants)



if __name__=="__main__":
    main()

