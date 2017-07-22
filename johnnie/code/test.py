import config, utils, snp_info

def test_make_toy():
  utils.make_toy(input_fn, toy_output_fn, toy_percentage)

def test_extract_exon():
  utils.extract_exon(input_fn, exon_output_fn)

def test_synthetic_dob():
  dobs = utils.synthetic_dob(2504)
  assert len(dobs) == 2504, "Error in utils.synthetic_dob()"

def test_shannon_and_simpson():
  snp_info.shannon_and_simpson(True)

if __name__ == '__main__':
  input_fn = '1000G_brca.hg38.vcf'
  toy_output_fn = input_fn[:-4]+'_toy.vcf'
  toy_percentage = 20
  exon_output_fn = input_fn[:-4]+'_exon.vcf'

  #test_make_toy()
  #test_extract_exon()
  #test_synthetic_dob()
  test_shannon_and_simpson()
