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

def read_snapout(snapout_name):
  """ Read SNAP output into dictionary.
      Each entry has a SNP being the key and all proxies and their LD info being the value.
  """
  snapout = csv.DictReader(open(config.data_path+snapout_name, 'rb'), delimiter='\t')
  ld_dict = {}
  for line in snapout:
    snp, prx, rsq, dpr = line['SNP'], line['Proxy'], line['RSquared'], line['DPrime']

    # Check if is error line
    if prx.startswith('WARNING') or prx.startswith('ERROR'):
      continue

    if snp in ld_dict.keys():
      lds = ld_dict[snp]
    else:
      lds = []
    lds.append({'proxy': prx, 'rsquared': rsq, 'dprime': dpr})
    ld_dict[snp] = lds

  return ld_dict


def write_snap_to_vcf(snapout_name, outvcf_name):
  """ Write a VCF file from SNAP output file. """
  # Read SNAP output as dicts
  snapout = csv.DictReader(open(config.data_path+snapout_name, 'rb'), delimiter='\t')
  for d in snapout:
    if d['Proxy'].startswith('WARNING') or d['Proxy'].startswith('ERROR'):
      pass

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
  # Don't need to: Convert from HG38 to HG18: 38 to 19, 19 to 18

  # Create query file
  write_query_rsids('1000G_brca.hg38to19to18_exon.vcf')

  # Manually query from SNAP by loading query file to their web app

  # Convert to VCF
  #write_snap_to_vcf('SNAPResults-chr13chr17-1K1.txt', 'test_writer.vcf')

