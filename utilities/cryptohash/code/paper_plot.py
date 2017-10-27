import snp_info
import snap_reader
import knapsack
import exprs
import utils
import config
import os
import pickle
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

input_vcf = '1000G_brca.hg38_exon.vcf'
snapout_name = 'SNAPResults-BRCA12exon-HapMap3r2_rsq08.txt'

def ranked_entropy():
  ranked = [r for r in range(1, len(ranked_snps)+1)]
  plt.plot(ranked, entropies, label='Entropy')
  snp_inds = [int(i)+1 for hap in hap_snp_inds for i in hap]
  print snp_inds
  plt.scatter(snp_inds, [entropies[i] for i in snp_inds], label='SNPs with ' + r'$r^2 > 0.8$')
  plt.title('SNP Entropies and Haplotypes')
  plt.legend(loc='best', shadow=True)
  plt.xlim([0, 550])
  plt.grid(True)
  plt.savefig(config.plot_path + 'exon_sorted_entropies_and_haps.pdf', bbox_inches='tight', pad_inches=0)

if __name__ == '__main__':
  ents_rankedsnps_path = config.data_path + 'ents_rankedsnps.pickle'
  if os.path.exists(ents_rankedsnps_path):
    print 'found', ents_rankedsnps_path
    entropies, ranked_snps = pickle.load(open(ents_rankedsnps_path, 'rb'))
  else:
    entropies, ranked_snps = snp_info.shannon_entropy(vcf.Reader(open(config.data_path + input_vcf, 'rb')))
    pickle.dump((entropies, ranked_snps), open(ents_rankedsnps_path, 'wb'))
  
  haps_path = config.data_path + 'big_haplotypes.pickle'
  if os.path.exists(haps_path):
    print 'found', haps_path
    hap_snp_inds =  pickle.load(open(haps_path, 'rb'))
  else:
    proxy_dict = snap_reader.read_snapout(snapout_name, input_vcf=input_vcf)
    hap_snp_inds = snap_reader.haplotype(proxy_dict)
    pickle.dump(hap_snp_inds, open(haps_path, 'wb'))

  """ 4.1 ranked entropy plot """
  ranked_entropy()




