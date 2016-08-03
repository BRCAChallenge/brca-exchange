import pysam
import overlapping_calculation_v2 as code

def main():
    v1 = ['13', '32890990', 'G', 'A']
    v2 = ['13', '32890990', 'G', 'A']
    print code.variant_pair_same(v1, v2)

def test_pysam_reference():
    f = pysam.FastaFile("references/chr13.fa")
    region = f.fetch(reference="chr13", start=20000000, end = 20000100)
    print region
    print type(region)
    print HAHA
    return region

def test_string_operations():
    s = test_pysam_reference()
    ref = "ta"
    alt = "gg"
    new_s = alt + s[len(ref):]
    print new_s
    new_s = s[:]


def test_exception():
    raise Exception("this just doesn't work")


def test_list_pop():
    a = generate_list()
    print a
    a.pop(2)
    print a



def test_append_lines_of_file():
    f = open("data/brca.clinvar.no_header.vcf")
    rows = []
    for line in f:
        rows.append(line)
        if len(rows) > 10:
            break
    print rows

def test_list_assignment():
    x, y = generate_list()[:2]
    print x, y
    print y

def test_variant_generation_and_swap():
    f = open("data/brca.clinvar.no_header.vcf")
    n = 0
    for line in f:
        n += 1
        if n == 1:
            chr1, pos1, ref1, alt1 = line.split("\t")[:2] + line.split("\t")[3:5]
        if n == 2:
            chr2, pos2, ref2, alt2 = line.split("\t")[:2] + line.split("\t")[3:5]
            break
    print chr1, pos1, ref1, alt1
    print chr2, pos2, ref2, alt2
    ref1, ref2 = ref2, ref1
    print chr1, pos1, ref1, alt1
    print chr2, pos2, ref2, alt2


def generate_list():
    return [1, 2, 3, 4, 5]



if __name__ == "__main__":
    main()
    
