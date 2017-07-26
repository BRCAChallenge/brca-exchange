import utils
import snp_info
import vcf
from matplotlib import pyplot as plt
import config
import numpy as np
import argparse

data_path = config.data_path
plot_path = config.plot_path

def init():
  parser = argparse.ArgumentParser(description="Specify parameters of experiment")
  parser.add_argument("-f", "--vcf_file", help="VCF file to use", default='1000G_brca.hg38_exon.vcf')
  #parser.add_argument("-e", "--experiments", help="File specifying experiments to run", default="def_exprs.txt")
  parser.add_argument("-s", "--snp_sampling", help="Run SNP sampling expr.", action='store_true')
  parser.add_argument("-a", "--af_plot", help="Run allele frequency expr.", action='store_true')
  parser.add_argument("-b", "--birth_type", help="Use date or year of birth", default='date')
  args = parser.parse_args()

  return args

def high_af_snps(high_af_fn):
  """  """
  pass

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
    # TODO: do more times and get mean
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

if __name__ == '__main__':
  args = init()

  file = args.vcf_file
  config.birth_type = args.birth_type

  if args.snp_sampling:
    snp_sampling(file)
  if args.af_plot:
    plot_af_per_snp(file)
