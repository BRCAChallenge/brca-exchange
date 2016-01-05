"""
this script add fields to a VCF file
including:
HGVS_G: HGVS genomic coordinate
HGVS_C: HGVS cDNA coordinate(s)
HGVS_P: HGVS protein coordinate(s)
gene_symbol
"""
# in vitae hgvs module is named hgvs, counsyl hgvs module is imported as pyhgvs
import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.variantmapper
import random
import repeat_merging
from pygr.seqdb import SequenceFileDB
import pyhgvs
import pyhgvs.utils as pyhgvs_utils


GENOME = SequenceFileDB('../data/hg19.fa')

def main():
    args = repeat_merging.arg_parse()
    add_HGVS_g(args.input, "../data/temp1")
    add_HGVS_c_counsyl("../data/temp1", "../data/temp2")
    add_HGVS_p_invitae("../data/temp2", args.output)

def add_HGVS_g(in_path, out_path):
    """
    counsyl pyhgvs module is used
    seems like invitae hgvs module doesn't have capacities to convert between
    hgvs genomic and vcf genomic format
    """
    f_in = open(in_path, "r")
    f_out = open(out_path, "w")
    line_num = 0
    for line in f_in:
        line_num += 1
        print line_num
        if line_num == 1:
            items = line.strip().split("\t")
            items.insert(3, "HGVS_genomic")
            new_line = "\t".join(items) + "\n"
        else:
            items = line.strip().split("\t")
            genome_coor = items[2].split(":")
            chrom = genome_coor[0]
            offset = int(genome_coor[1])
            ref = genome_coor[2].split(">")[0]
            alt = genome_coor[2].split(">")[1]
            hgvs_name = pyhgvs.variant_to_hgvs_name(chrom, offset, ref, alt, GENOME, None)
            hgvs_g = "NC_0000" + chrom[-2:] + ".10:g." + hgvs_name.format_genome()
            items.insert(3, hgvs_g)
            new_line = "\t".join(items) + "\n"
        f_out.write(new_line)


def add_gene_symbol_and_transcript(in_path, out_path):
    f_in = open(in_path, "r")
    f_out = open(out_path, "w")
    line_num = 0
    for line in f_in:
        line_num += 1
        print line_num
        if line_num == 1:
            f_out.write(line)
        else:
            items = line.strip().split("\t")
            chrm = items[2][:5]
            if chrm == "chr13":
                gene_symbol = "BRCA2"
                transcript = "NM_000059.3"
            elif chrm == "chr17":
                gene_symbol = "BRCA1"
                transcript = "NM_007294.3"
            else:
                raise Exception("no genome coordinate or malformed genome coordinate")

            # add gene symbol, BRCA1 for chromsome 17, BRCA2 for chromosome 13
            if items[1] == "-":
                items[1] = gene_symbol
            else:
                if items[1] == gene_symbol:
                    pass
                else:
                    raise Exception("mismatched gene symbol")

            # add cDNA transcript id
            # NM_000059.3 for chromosome 13, NM_007294.3 for chromosome 17
            if items[3] == "-":
                items[3] = transcript
            else:
                if items[3] == transcript:
                    pass
                else:
                    raise Exception("mismatched transcript")

            new_line = "\t".join(items) + "\n"
            f_out.write(new_line)


def add_info(line):
    items = line.strip().split("\t")
    current_info = items[-1]
    current_info += ";FREQ=" + str(random.random())
    current_info += ";BIC_P=?;BIC_N=?"
    current_info += ";MUTTYPE=" + random.choice(["missense", "nonsense", "silent"])
    current_info += ";IARC=" + random.choice(["group 1", "group 2a", "group 2b", "group 3", "group 4"])
    current_info += ";DBSource=" + random.choice(["Clinvar", "LOVD", "BIC", "ExAC", "1000Genomes"])
    items[-1] = current_info
    new_line = "\t".join(items) + "\n"
    return new_line


def add_HGVS_c_counsyl(in_path, out_path):
    # counsyl pyhgvs setups
    with open('../data/BRCA12.refGene.txt') as infile:
        transcripts_counsyl = pyhgvs_utils.read_transcripts(infile)
    def get_transcript(name):
        return transcripts_counsyl.get(name)

    f_in = open(in_path, "r")
    f_out = open(out_path, "w")
    line_num = 0
    unmatching_cases = 0
    for line in f_in:
        line_num += 1
        print line_num
        if line_num == 1:
            f_out.write(line)
            continue
        items = line.strip().split("\t")
        genome_coors = items[2].split(":")
        chrom = genome_coors[0]
        offset = int(genome_coors[1])
        ref = genome_coors[2].split(">")[0]
        alt = genome_coors[2].split(">")[1]
        transcript_id = items[4]
        transcript = get_transcript(transcript_id)
        hgvs_name_with_transcript = pyhgvs.variant_to_hgvs_name(
            chrom, offset, ref, alt, GENOME, transcript)
        hgvs_c = hgvs_name_with_transcript.format(use_prefix=False,
                                                  use_gene=False)
        if items[5] == "-":
            items[5] = hgvs_c
        elif items[5] == hgvs_c:
            pass
        else:
            unmatching_cases += 1
        new_line = "\t".join(items) + "\n"
        f_out.write(new_line)
    print "unmatching cases: ", unmatching_cases

def add_HGVS_c_invitae(in_path, out_path):
    f_in = open(in_path, "r")
    f_out = open(out_path, "w")
    # in vitae hgvs setups
    hp = hgvs.parser.Parser()
    hdp = hgvs.dataproviders.uta.connect()
    evm = hgvs.variantmapper.EasyVariantMapper(hdp,
        primary_assembly = 'GRCh37', alt_aln_method = 'splign')
    line_num = 0
    unmatching_cases = 0
    not_implemented_cases = 0
    for line in f_in:
        line_num += 1
        print line_num
        if line_num == 1:
            f_out.write(line)
            continue
        items = line.strip().split("\t")
        HGVS_genomic = items[3]
        try:
            hgvs_g = hp.parse_hgvs_variant(HGVS_genomic)
            transcript_id = items[4]
            try:
                hgvs_c = str(evm.g_to_c(hgvs_g, transcript_id)).split(":")[1]
            except:
                hgvs_c = "not implemented"
                not_implemented_cases += 1
        except hgvs.exceptions.HGVSParseError:
            not_implemented_cases += 1
            hgvs_c = "not implemented"

        if items[5] == "-":
            items[5] = hgvs_c
        elif items[5] == hgvs_c:
            pass
        elif hgvs_c != "not implemented":
            unmatching_cases += 1
        new_line = "\t".join(items) + "\n"
        f_out.write(new_line)
    print "unmatching cases: ", unmatching_cases
    print "not implemented cases: ", not_implemented_cases

def compare_counsyl_and_invitae(counsyl_file, invitae_file):
    f_c = open(counsyl_file, "r")
    f_i = open(invitae_file, "r")
    clist = []
    ilist = []
    unmatching_num = 0
    unmatching_implemented_num = 0
    for line in f_c:
        clist.append(line.strip())
    for line in f_i:
        ilist.append(line.strip())
    for i in range(len(clist)):
        if clist[i] != ilist[i]:
            unmatching_num += 1
            if ilist[i] != "not implemented":
                unmatching_implemented_num += 1
                print clist[i]
                print ilist[i]
                print ""
    print "unmatching cases between counsyl and invitae: ", unmatching_num
    print "unmatching implemented cases between counsyl and invitae: ", unmatching_implemented_num

def add_HGVS_p_invitae(in_path, out_path):
    # in vitae hgvs setups
    hp = hgvs.parser.Parser()
    hdp = hgvs.dataproviders.uta.connect()
    evm = hgvs.variantmapper.EasyVariantMapper(hdp,
        primary_assembly = 'GRCh37', alt_aln_method = 'splign')

    f_in = open(in_path, "r")
    f_out = open(out_path, "w")
    line_num = 0
    unmatching_cases = 0
    notimplemented_cases = 0
    for line in f_in:
        line_num += 1
        print line_num
        if line_num == 1:
            f_out.write(line)
            continue
        items = line.strip().split("\t")
        HGVS_cDNA_string = items[4] + ":" + items[5]
        try:
            hgvs_c = hp.parse_hgvs_variant(HGVS_cDNA_string)
            hgvs_p = str(evm.c_to_p(hgvs_c)).split(":")[1]
        except:
            notimplemented_cases += 1
            hgvs_p = "-"

        if items[6] == "-":
            items[6] = hgvs_p
        elif items[6] == hgvs_p:
            pass
        else:
            unmatching_cases += 1
        new_line = "\t".join(items) + "\n"
        f_out.write(new_line)
    print "unmatching cases: ", unmatching_cases
    print "notimplemented cases: ", notimplemented_cases



if __name__ == "__main__":
    main()






