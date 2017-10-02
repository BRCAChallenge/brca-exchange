import vcf
import snp_info
import snap_reader
import config
import argparse
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import utils
import pickle
import os

input_vcf = '1000G_brca.hg38_exon.vcf'
ents_rankedsnps_path = config.data_path + 'ents_rankedsnps.pickle'
if os.path.exists(ents_rankedsnps_path):
  print 'found', ents_rankedsnps_path
  entropies, ranked_snps = pickle.load(open(ents_rankedsnps_path, 'rb'))
else:
  entropies, ranked_snps = snp_info.shannon_entropy(vcf.Reader(open(config.data_path + input_vcf, 'rb')))
  pickle.dump((entropies, ranked_snps), open(ents_rankedsnps_path, 'wb'))

snapout = 'SNAPResults-BRCA12exon-HapMap3r2_rsq08.txt'
proxy_dict_path = config.data_path + 'proxy_dict.pickle'
if os.path.exists(proxy_dict_path):
  print 'found', proxy_dict_path
  proxy_dict = pickle.load(open(proxy_dict_path, 'rb'))
else:
  proxy_dict = snap_reader.read_snapout(snapout, input_vcf=input_vcf)
  pickle.dump(proxy_dict, open(proxy_dict_path, 'wb'))

hap_inds_path = config.data_path + 'big_haplotypes.pickle'
if os.path.exists(hap_inds_path):
  print 'found', hap_inds_path
  hap_inds = pickle.load(open(hap_inds_path, 'rb'))
else:
  big_haplotypes = [hap for hap in snap_reader.haplotype(proxy_dict) if len(hap) > 1]
  hap_inds = [sorted([ranked_snps.index(snp) for snp in hap]) for hap in big_haplotypes]
  pickle.dump(hap_inds, open(hap_inds_path, 'wb'))


def init():
  parser = argparse.ArgumentParser(description='Specify parameters of experiment')
  parser.add_argument('-e', '--experiment', help='Experiment to run:\n1. Vary C on x-axis, run W times\n2. Vary W on x-axis, run C times', type=int, default=1)
  parser.add_argument('-w', '--maxweight', help='Maximum weight of selection', type=int, default=1)
  parser.add_argument('-c', '--count', help='Count of items to select', type=int, default=25)
  args = parser.parse_args()
  return args

def solve_knapsack(W, C):
  weights = []
  values = []
  N = len(ranked_snps)
  for snp_i in range(N):
    snp = ranked_snps[snp_i]
    if snp in proxy_dict.keys():
      proxies = proxy_dict[snp]
    else:
      proxies = []
    weights.append(weight(snp_i, N, proxies))
    values.append(value(snp_i))

  # Normalize weights to have sum 100
  w_sum = sum(weights)
  weights = [w*100/w_sum for w in weights]

  for hap in hap_inds:
    num = len(hap)/2
    for ind in hap:
      if num == 0:
        continue
      weights[ind] = 0
      num -= 1

  #for i in range(30):
  #  weights[i] = 0

  return knapsack(weights, values, W, N, C)

def knapsack(ws, vs, W, N, C):
  """ Dynamic Programming solution for 0/1 Knapsack problem.
      Returns the maximum value that can be put in a knapsack of capacity W
      and the items that give that value.
  """
  M = [[[(0, []) for k in range(C+1)] for x in range(W+1)] for x in range(N)]
  for i in range(N):
    for j in range(W+1):
      for k in range(C+1):
        if ws[i] > j or k < 1:
          M[i][j][k] = M[i-1][j][k]
        else:
          prev_val, prev_its = M[i-1][int(j-ws[i])][k-1]
          leav_val, leav_its = M[i-1][j][k]
          take_val = vs[i] + prev_val
          take_its = prev_its + [ranked_snps[i]]

          M[i][j][k] = max((take_val, take_its), (leav_val, leav_its))

  return M[N-1][W][C]

def weight(snp_i, snp_num, proxies):
  """ Proportion of maximum possible """
  return sum(r2(snp_i, p) for p in proxies) / snp_num

def r2(snp_i, proxy):
  """ Return correlation value. """
  return float(proxy['rsquared'])

def value(snp_i):
  """ Return entropy of SNP. """
  return entropies[snp_i]

def expr_vary_W(W, C):
  idabs = []
  config.birth_type = 'year'
  dobs = utils.synthetic_dob(2504)
  dobs = None
  for w in range(1, W+1):
    value, knap_snps = solve_knapsack(w, C)
    #idab, _ = snp_info.identifiability(input_vcf, True, dobs=dobs, snps=knap_snps)
    idab, _ = snp_info.identifiability(input_vcf, False, dobs=dobs, snps=knap_snps)
    print 'accuracy at', w, 'th round:', idab
    idabs.append(idab)
  return idabs

def expr_vary_C(W, C, crange):
  idabs = []
  for c in crange:
    value, knap_snps = solve_knapsack(W, c)
    idab, _ = snp_info.identifiability(input_vcf, False, snps=knap_snps)
    print 'accuracy at', c, 'th round:', idab
    idabs.append(idab)
  return idabs

def plot_expr_results(idabs, var_WC, var, title=None, xlabel=None, xticks=None):
  plt.plot(range(var_WC), idabs)
  if title:
    plt.title(title)
  if xlabel:
    plt.xlabel(xlabel)
  if xticks:
    plt.xticks(xticks)
  pickle.dump(idabs, open(config.data_path + 'knapsack_vary_' + var + '_' + str(var_WC) + '.pickle', 'wb'))
  plt.savefig(config.plot_path + 'knapsack_vary_' + var + '_' + str(var_WC) + '.pdf')

def plot_mult_expr_results(idabs_lst, xvar_len, ovar_range, figname='', title='', xlabel='', xticks=None):
  for i, idabs in enumerate(idabs_lst):
    plt.plot(range(len(idabs)), idabs, label='W = '+str(ovar_range[i]))
  if title:
    plt.title(title)
  if xlabel:
    plt.xlabel(xlabel)
  if xticks:
    plt.xticks(range(xvar_len), xticks)
  plt.legend(loc='best', shadow=True)
  plt.savefig(config.plot_path + figname)

if __name__ == '__main__':
  args = init()

  W = args.maxweight
  C = args.count

  if args.experiment == 1:
    idabs_lst_path = config.data_path+'knapsack_x_C'+str(C)+'_W'+str(W)+'.pickle'
    figname = 'knapsack_x_C'+str(C)+'_W'+str(W) + '.pdf'
    title = 'Identifiability with Varying Number of SNPs'
    xlabel = 'Number of SNPs'
  elif args.experiment == 2:
    idabs_lst_path = config.data_path+'knapsack_x_W'+str(W)+'_C'+str(C)+'.pickle'
    figname = 'knapsack_x_W'+str(W)+'_C'+str(C) + '.pdf'
    title = 'Identifiability with Varying Max. Weight'
    xlabel = 'Maximum Weight'

  if W >= 10:
    wsteps = range(0, W, W/10)
  else:
    wsteps = range(W)
  if C >= 10:
    crange = range(1, C+1, (C+1)/10)
  else:
    crange = range(1, C+1)

  if os.path.exists(idabs_lst_path):
    idabs_lst = pickle.load(open(idabs_lst_path, 'rb'))
  else:
    idabs_lst = []
    for w in wsteps:
      idabs = expr_vary_C(w, C, crange)
      idabs_lst.append(idabs)
      print 'w', w, 'C', C, 'num of idabs', len(idabs)
      print
    pickle.dump(idabs_lst, open(idabs_lst_path, 'wb'))

  plot_mult_expr_results(idabs_lst, len(crange), wsteps, figname=figname, title=title, xlabel=xlabel, xticks=[str(r) for r in crange])

  #print idabs
  #plot_expr_results(idabs, W, 'W', title='Identifiability with 30 SNPs', xlabel='Max. Weight')
  #plot_expr_results(idabs, C, 'C', title='Identifiability with Max. Weight '+str(w), xlabel='Number of SNPs')
