import utilities
import snp_info
import vcf
from matplotlib import pyplot as plt
import config
import numpy as np

data_path = config.data_path
plot_path = config.plot_path

def snp_sampling(exon_fn):
  ''' Vary the size of subsets of SNPs (10%, 20%, ..., 100%), 
      and show how uniqueness is affected. '''
  subset_sizes = [10 * (i+1) for i in range(10)]
  sub_fn = exon_fn[:-4] + '.subset.vcf'
  uniqs_dob = []
  uniqs = []
  dobs = utilities.synthetic_dob(2504)
  for size in subset_sizes:
    print 'Using random', size, '% of SNPs'
    # TODO: do more times and get mean
    utilities.make_toy(data_path+exon_fn, data_path+sub_fn, size)
    uniq_dob = snp_info.uniqueness(data_path+sub_fn, use_dob=True, dobs=dobs)
    uniq = snp_info.uniqueness(data_path+sub_fn, use_dob=False, dobs=dobs)
    uniqs_dob.append(uniq_dob)
    uniqs.append(uniq)
    print

  plt.plot(subset_sizes, uniqs, 'b', label='w/o DOB')
  plt.plot(subset_sizes, uniqs_dob, 'r', label='w/  DOB')
  ax = plt.subplot()
  d = (subset_sizes[1] - subset_sizes[0]) * 0.1
  xlim = [subset_sizes[0] - d, subset_sizes[-1] + d]
  ylim = [-0.05, 1.05]
  ax.set_xlim(xlim)
  ax.set_ylim(ylim)
  plt.legend(loc='best', shadow=True)
  plt.title('Identifiability using random X % of SNPs')
  plt.savefig(plot_path+'snp_sampling.pdf')

def plot_af_per_snp(input_fn):
  vcf_reader = vcf.Reader(open(data_path+input_fn, 'rb'))
  afs = []
  max_a_num = 0
  for snp in vcf_reader:
    af = utilities.allele_frequencies(snp)
    afs.append(af)
    if len(af) > max_a_num:
      max_a_num = len(af)

  afs = utilities.fill_zeros(afs, max_a_num)
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
  plt.savefig(plot_path+'allele_frequency_per_SNP.pdf')

if __name__ == '__main__':
  #snp_sampling('1000G_brca.hg38_exon.vcf')
  plot_af_per_snp('1000G_brca.hg38_exon.vcf')
