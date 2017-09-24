import utils
import snp_info
import vcf
from matplotlib import pyplot as plt
import config
import numpy as np
import argparse
import hashlib

data_path = config.data_path
plot_path = config.plot_path

def init():
  parser = argparse.ArgumentParser(description="Specify parameters of experiment")
  parser.add_argument("-i", "--input", help="Input VCF file.", default='1000G_brca.hg38_exon.vcf')
  parser.add_argument("-o", "--output", help="Output VCF file.", default='output.txt')
  parser.add_argument("-s", "--snp_sampling", help="Run SNP sampling expr.", action='store_true')
  parser.add_argument("-a", "--af_plot", help="Run allele frequency expr.", action='store_true')
  parser.add_argument("-t", "--af_threshold", help="Run allele freq. threshold expr.", default=None)
  parser.add_argument("-c", "--crypto_hash", help="Generate cryptographic hash.", action='store_true')
  parser.add_argument("-b", "--birth_type", help="Use date or year of birth", default='')
  args = parser.parse_args()

  return args

def high_af_snps(high_af_fn):
  """  """
  pass

def expr_vary_af_treshold(input_fn, type='a'):
  if type == 'a': # 0, 1, ..., 10
    thresholds = [r for r in range(0, 11)]
  if type == 'b': # 0, 0.1, ..., 1
    thresholds = [0.1 * r for r in range(0, 11)]
  plot_vary_af_threshold(input_fn, thresholds, type)

def plot_vary_af_threshold(input_fn, thresholds, type):
  """ Vary allele frequency threshold and plot change in identifiability. """
  dobs = utils.synthetic_dob(2504)

  # Find identifiability of SNPs selected using each threshold.
  idabs = []
  snp_nums = []
  for t in thresholds:
    af_thresh_fn = input_fn[:-4] + '_af' + str(t) + '.vcf'
    utils.extract_high_af(input_fn, af_thresh_fn, t)
    idab, snp_num = snp_info.identifiability(af_thresh_fn, use_dob=True, dobs=dobs)
    idabs.append(idab)
    snp_nums.append(snp_num)

  # Plot identifiability change.
  plt.plot(thresholds, idabs)
  # In xticks show number of SNPs obtained with each threshold
  plt.xticks(thresholds, [t+'\n'+str(snp_nums[i]) for i, t in enumerate(thresholds)])
  plt.title('Identifiability using SNPs with each allele frequency >= t%')
  plt.savefig(plot_path+'thresholds'+'_'+type+'_birth'+config.birth_type+'.pdf')

def snp_sampling(input_fn):
  """ Vary the size of subsets of SNPs (10%, 20%, ..., 100%), 
      and show how identifiability is affected. """
  subset_sizes = [10 * (i+1) for i in range(10)]
  sub_fn = input_fn[:-4] + '.subset.vcf'
  uniqs_dob = []
  uniqs = []
  dobs = utils.synthetic_dob(2504)
  for size in subset_sizes:
    print 'Using random', size, '% of SNPs'
    utils.make_toy(data_path+input_fn, data_path+sub_fn, size)
    uniq_dob, snp_num = snp_info.identifiability(data_path+sub_fn, use_dob=True, dobs=dobs)
    uniq, snp_num = snp_info.identifiability(data_path+sub_fn, use_dob=False, dobs=dobs)
    uniqs_dob.append(uniq_dob)
    uniqs.append(uniq)
    print

  if config.birth_type == 'date':
    dob = 'DOB'
  elif config.birth_type == 'year':
    dob = 'YOB'
  plt.plot(subset_sizes, uniqs, 'b', label='w/o '+dob)
  plt.plot(subset_sizes, uniqs_dob, 'r', label='w/ '+dob)
  ax = plt.subplot()
  d = (subset_sizes[1] - subset_sizes[0]) * 0.1
  xlim = [subset_sizes[0] - d, subset_sizes[-1] + d]
  ylim = [-0.05, 1.05]
  ax.set_xlim(xlim)
  ax.set_ylim(ylim)
  plt.legend(loc='best', shadow=True)
  plt.title('Identifiability using random X % of '+ str(snp_num) +' SNPs')
  plt.savefig(plot_path+'snp_sampling_'+input_fn[:-4]+'_'+dob+'.pdf')
  plt.close()

def plot_af_per_snp(input_fn):
  """ Generate allele frequency plot, sorted by. """
  vcf_reader = vcf.Reader(open(data_path+input_fn, 'rb'))
  afs = []
  max_a_num = 0
  for snp in vcf_reader:
    af = utils.allele_frequencies(snp)
    afs.append(af)
    if len(af) > max_a_num:
      max_a_num = len(af)

  afs = utils.fill_zeros(afs, max_a_num)
  afs = sorted(afs, key=lambda i: i[0], reverse=True)
  afs = np.array(afs).T
  ind = np.arange(len(afs[0]))
  prevtop = np.zeros(len(afs[0]))
  colors = ['b', 'r', 'g', 'b', 'm', 'c']
  for i, af in enumerate(afs):
    af = np.array(af) + np.array(prevtop)
    plt.plot(ind, af, color=colors[i], alpha=0.5, label=str(i))
    prevtop = af
  ax = plt.subplot()
  ax.set_ylim([-0.05, 1.05])
  plt.title('Allele frequency sorted by highest frequency.')
  plt.legend(loc='best', shadow=True)
  plt.savefig(plot_path+'af_per_SNP_'+input_fn[:-4]+'.pdf')
  plt.close()

def hash_func(strings):
  """ Function to compute hashes, using hexdigest of SHA256. """
  return [hashlib.sha256(s).hexdigest() for s in strings]

def make_hash(input_fn, output_fn, use_dob=True, dobs=None):
  """ Generate hash code for individuals in input file. """
  vcf_reader = vcf.Reader(open(data_path+input_fn, 'rb'))
  m = len(vcf_reader.samples)
  n = sum(1 for _ in vcf_reader)
  if dobs == None:
    dobs = utils.synthetic_dob(m)
  genotypes = snp_info.genotypes(input_fn, m, n, use_dob, dobs)

  hashes = hash_func(genotypes)
  with open(data_path+output_fn, 'wb') as outfile:
    for hash in hashes:
      outfile.write(hash)
      outfile.write('\n')

if __name__ == '__main__':
  args = init()

  infile = args.input
  outfile = args.output
  config.birth_type = args.birth_type

  if args.snp_sampling:
    snp_sampling(infile)
  if args.af_plot:
    plot_af_per_snp(infile)
  if args.crypto_hash:
    make_hash(infile, outfile)
  if args.af_threshold:
    expr_vary_af_treshold(infile, args.af_threshold)
