import vcf
import csv
import config

def write_query_rsids(infile_name, outfile_name=None):
  """ Write RSIDs of all SNPs in infile (VCF) to outfile, one on each line. """
  assert infile_name[-4:] == '.vcf', 'infile should be a VCF file'
  vcf_reader = vcf.Reader(open(config.data_path+infile_name, 'rb'))
  if not outfile_name:
    outfile_name = infile_name[:-4] + '_ids.txt'
  with open(config.data_path+outfile_name, 'wb') as outfile:
    for snp in vcf_reader:
      outfile.write(snp.ID+'\n')

def write_snap_to_vcf(snapout_name, outvcf_name):
  """ Write a VCF file from SNAP output file. """
  # Read SNAP output as dicts
  snapout = csv.DictReader(open(config.data_path+snapout_name, 'rb'), delimiter='\t')
  for d in snapout:
    if d['Proxy'].startswith('WARNING') or d['Proxy'].startswith('ERROR'):
      pass

  # Read vcf as dicts
  #vcf_reader = vcf.Reader(open(config.data_path+'1000G_brca.hg38_exon.vcf', 'rb'))
  #myout = []
  #for snp in vcf_reader:
  #  d = {}
  #  d['Coordinate_HG38'] = snp.POS
  #  d['Chromosome'] = snp.CHROM
  #  d['SNP'] = snp.ID
  #  myout.append(d)

  # Create VCF file from SNAP
  vcf_reader = vcf.Reader(open(config.data_path+'1000G_brca.hg38_exon.vcf', 'rb'))
  vcf_writer = vcf.Writer(open(config.data_path+outvcf_name, 'wb'), vcf_reader)
  for d in snapout:
    print d
    try:
      record = vcf.model._Record(CHROM=d['Chromosome'].split('chr')[1], POS=int(d['Coordinate_HG18']), ID=d['Proxy'], REF=[], ALT=[], QUAL=[], FILTER=[], INFO=[], FORMAT=[], sample_indexes=[])
      vcf_writer.write_record(record)
    except AttributeError:
      pass
  vcf_writer.close()


if __name__ == '__main__':
  # Create query file
  #write_query_rsids('1000G_brca.hg38_exon.vcf')

  # Manually query from SNAP

  # Convert to VCF
  write_snap_to_vcf('SNAPResults-chr13chr17-1K1.txt', 'test_writer.vcf')
