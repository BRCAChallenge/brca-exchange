import matplotlib
matplotlib.use('Agg')
import utils
import snp_info
import vcf
from matplotlib import pyplot as plt
import config
import numpy as np
import argparse
import hashlib
import os
import pickle

data_path = config.data_path
plot_path = config.plot_path

def init():
  parser = argparse.ArgumentParser(description='Specify parameters of experiment')
  parser.add_argument('-i', '--input', help='Input VCF file.', default='1000G_brca.hg38_exon.vcf')
  parser.add_argument('-o', '--output', help='Output VCF file.', default='output.txt')
  parser.add_argument('-s', '--snp_sampling', help='Run SNP sampling expr.', action='store_true')
  parser.add_argument('-a', '--af_plot', help='Run allele frequency expr.', action='store_true')
  parser.add_argument('-t', '--af_threshold', help='Run AF threshold expr.: 1-3', type=int, default=0)
  parser.add_argument('-d', '--independence', help='Run independence expr.', type=int, default=0)
  parser.add_argument('-e', '--entropy', help='Run top k entropy expr.: 1-3', type=int, default=0)
  parser.add_argument('-f', '--entind', help='Run entropy-independence expr.', type=int, default=0)
  parser.add_argument('-c', '--crypto_hash', help='Generate cryptographic hash.', action='store_true')
  parser.add_argument('-b', '--birth_type', help='Use date or year of birth', default='')
  parser.add_argument('-p', '--snap_output', help='SNAP output .txt file', default='SNAPResults-BRCA12exon-HapMap3r2_rsq08.txt')
  args = parser.parse_args()

  return args

def expr_top_k_independent(input_fn):
  """  """
  ks = range(1, 51)
  births = ['date', 'year', '']
  plot_data = []
  for birth in births:
    print birth
    data = get_data_top_k_independent(input_fn, ks, birth, 1)
    plot_data.append(data)
  plot_top_k_entropy2(ks, plot_data, births, 1)

def plot_top_k_entropy2(ks, plot_data, births, expr):
  colors = ['r', 'b', 'g', 'm', 'c']
  ymax, ymin = 0.0, 1.0
  for i, (idabs, snp_nums) in enumerate(plot_data):
    birth_type = births[i] if births[i] else 'no birth'
    plt.plot(ks, idabs, colors[i], label=birth_type, alpha=1)
    ymax = max(ymax, max(idabs))
    ymin = min(ymin, min(idabs))
  plt.legend(loc='best', shadow=True)
  plt.title('Identifiability of Top k Independent SNPs Ranked by Entropy')
  plt.xlabel('Number k of SNPs')
  xmax, xmin = max(ks), min(ks)
  xrng = xmax - xmin
  yrng = ymax - ymin
  plt.xlim([xmin - 0.05 * xrng, xmax + 0.05 * xrng])
  plt.ylim([ymin - 0.05 * yrng, ymax + 0.05 * yrng])
  plt.xlim([xmin, xmax])
  plt.ylim([0, 1.02])
  plt.grid(True)
  plt.savefig(plot_path + 'top_k_entropy_independence.pdf', bbox_inches='tight', pad_inches=0)

def get_data_top_k_independent(input_fn, ks, birth, expr):
  if birth:
    config.birth_type = birth
    dobs = utils.synthetic_dob(2504)
  else:
    dobs = []
  return get_idabs_top_k_independent(input_fn, ks, birth, dobs, expr)

def get_idabs_top_k_independent(input_fn, ks, birth, dobs, expr):
  birth_type = birth if birth else 'nobirth'
  expr_data_path = data_path + 'top_k_indep_'+birth_type+'.pickle'
  if os.path.exists(expr_data_path):
    idabs, snp_nums = pickle.load(open(expr_data_path, 'rb'))
  else:
    idabs = []
    snp_nums = []
    rsq = 8 # 0.8
    proxy_fn = 'SNAPResults-BRCA12exon-HapMap3r2_rsq'+ str(rsq).zfill(2) +'.txt'
    for k in ks:
      selection = utils.extract_low_redund_2(input_fn, proxy_fn, k=k)
      idab, snp_num = snp_info.identifiability(input_fn, use_dob=birth, dobs=dobs, snps=selection)
      idabs.append(idab)
      snp_nums.append(snp_num)
    pickle.dump((idabs, snp_nums), open(expr_data_path, 'wb'))
  return idabs, snp_nums

def expr_vary_rsq(input_fn):
  """ Run experiment to compare identifiabilities of SNPs
      reduced from using varying r2 thresholds.
  """
  rsqs = range(1, 11)
  births = ['date', 'year', '']
  plot_data = []
  for birth in births:
    print birth
    data = get_data_vary_rsq(input_fn, rsqs, birth, 1)
    plot_data.append(data)
  plot_vary_rsq(rsqs, plot_data, births, 1)

def plot_vary_rsq(rsqs, plot_data, births, expr):
  colors = ['r', 'b', 'g', 'm', 'c']
  ymax, ymin = 0.0, 1.0
  for i, (idabs, snp_nums) in enumerate(plot_data):
    birth_type = births[i] if births[i] else 'no birth'
    plt.plot(rsqs, idabs, colors[i], label=birth_type, alpha=1)
    ymax = max(ymax, max(idabs))
    ymin = min(ymin, min(idabs))
  plt.legend(loc='best', shadow=True)
  plt.title('Identifiability of SNPs After Removing Correlated SNPs with Varying r^2')
  plt.xticks(rsqs, [str(0.1*rsq)+'\n'+str(snp_nums[i]) for i, rsq in enumerate(rsqs)])
  plt.xlabel('r^2')
  xmax, xmin = max(rsqs), min(rsqs)
  xrng = xmax - xmin
  yrng = ymax - ymin
  plt.xlim([xmin - 0.05 * xrng, xmax + 0.05 * xrng])
  plt.ylim([ymin - 0.05 * yrng, ymax + 0.05 * yrng])
  plt.savefig(plot_path + 'vary_rsq.pdf')  

def get_data_vary_rsq(input_fn, rsqs, birth, expr):
  if birth:
    config.birth_type = birth
    dobs = utils.synthetic_dob(2504)
  else:
    dobs = []
  print len(dobs)
  return get_idabs_vary_rsq(input_fn, rsqs, birth, dobs, expr)

def get_idabs_vary_rsq(input_fn, rsqs, birth, dobs, expr):
  birth_type = birth if birth else 'nobirth'
  expr_data_path = data_path + 'vary_rsq_data_'+birth_type+'.pickle'
  if os.path.exists(expr_data_path):
    idabs, snp_nums = pickle.load(open(expr_data_path, 'rb'))
  else:
    idabs = []
    snp_nums = []
    selections = []
    for rsq in rsqs:
      proxy_fn = 'SNAPResults-BRCA12exon-HapMap3r2_rsq'+ str(rsq).zfill(2) +'.txt'
      selection = utils.extract_low_redund_2(input_fn, proxy_fn)
      idab, snp_num = snp_info.identifiability(input_fn, use_dob=birth, dobs=dobs, snps=selection)
      idabs.append(idab)
      snp_nums.append(snp_num)
      selections.append(selection)
    pickle.dump((idabs, snp_nums), open(expr_data_path, 'wb'))
  return idabs, snp_nums
  
def expr_top_k_entropy(input_fn):
  """ Find identifiability of top k SNPs ranked with entropy, varying k. """
  ks = range(1, 51)
  births = ['date', 'year', '']
  plot_data = []
  for birth in births:
    print birth
    data = get_data_top_k_entropy(input_fn, ks, birth, 1)
    plot_data.append(data)
  plot_top_k_entropy(ks, plot_data, births, 1)

def plot_top_k_entropy(ks, plot_data, births, expr):
  colors = ['r', 'b', 'g', 'm', 'c']
  ymax, ymin = 0.0, 1.0
  for i, idabs in enumerate(plot_data):
    birth_type = births[i] if births[i] else 'no birth'
    plt.plot(ks, idabs, colors[i], label=birth_type, alpha=1)
    ymax = max(ymax, max(idabs))
    ymin = min(ymin, min(idabs))
  plt.legend(loc='best', shadow=True)
  plt.title('Identifiability of Top k SNPs Ranked by Entropy')
  plt.xlabel('Number k of SNPs')
  xmax, xmin = max(ks), min(ks)
  xrng = xmax - xmin
  yrng = ymax - ymin
  plt.xlim([xmin - 0.05 * xrng, xmax + 0.05 * xrng])
  plt.ylim([ymin - 0.05 * yrng, ymax + 0.05 * yrng])
  plt.xlim([xmin, xmax])
  plt.ylim([0, 1.02])
  plt.grid(True)
  plt.savefig(plot_path + 'top_k_entropy.pdf', bbox_inches='tight', pad_inches=0)

def get_data_top_k_entropy(input_fn, ks, birth, expr):
  if birth:
    config.birth_type = birth
    dobs = utils.synthetic_dob(2504)
  else:
    dobs = []
  print len(dobs)
  return get_idabs_top_k_entropy(input_fn, ks, birth, dobs, expr)

def get_idabs_top_k_entropy(input_fn, ks, birth, dobs, expr):
  birth_type = birth if birth else 'nobirth'
  #expr_path = 'top_k_entropy'
  expr_data_path = data_path + 'top_k_entropy_data_'+birth_type+'.pickle'
  if os.path.exists(expr_data_path):
    idabs = pickle.load(open(expr_data_path, 'rb'))
  else:
    vcf_reader = vcf.Reader(open(data_path + input_fn, 'rb'))
    entropies, snpids = snp_info.shannon_entropy(vcf_reader)
    #print entropies
    idabs = []
    for k in ks:
      topk_entropies, topk_snps = entropies[:k], snpids[:k]
      idab, _ = snp_info.identifiability(input_fn, use_dob=birth, dobs=dobs, snps=topk_snps)
      idabs.append(idab)
    pickle.dump(idabs, open(expr_data_path, 'wb'))
  return idabs

def expr_vary_af_treshold(input_fn, expr):
  thresholds = [r for r in range(0, 11)]
  births = ['date', 'year', '']
  plot_data = []
  for birth in births:
    data = get_data_vary_af_threshold(input_fn, thresholds, birth, expr)
    plot_data.append(data)
  plot_vary_af_threshold(thresholds, plot_data, births, expr)

def plot_vary_af_threshold(thresholds, plot_data, births, expr):
  # Plot identifiability change.
  colors = ['r', 'b', 'g', 'm', 'c']
  ymax, ymin = 0.0, 1.0
  for i, (idabs, snp_nums) in enumerate(plot_data):
    birth_type = births[i] if births[i] else 'no birth'
    plt.plot(thresholds, idabs, colors[i], label=birth_type, alpha=1)
    ymax = max(ymax, max(idabs))
    ymin = min(ymin, min(idabs))

  # In xticks show number of SNPs obtained with each threshold
  plt.xticks(thresholds, [str(t)+'%\n'+str(snp_nums[i]) for i, t in enumerate(thresholds)])
  plt.legend(loc='best', shadow=True)

  title = 'Identifiability of SNPs Filtered with '
  af_ttl = 'min AF Threshold x% '
  snap_ttl = 'r2 Threshold 0.8 '
  then_ttl = '\nand then '
  if expr == 1:
    title = (title + af_ttl).strip()
  elif expr == 2:
    title = (title + af_ttl + then_ttl + snap_ttl).strip()
  elif expr == 3:
    title = (title + snap_ttl + then_ttl + af_ttl).strip()
  plt.title(title)

  plt.xlabel('AF Threshold % and SNP Number')
  xmax, xmin = max(thresholds), min(thresholds)
  xrng = xmax - xmin
  ymin = 0.45
  yrng = ymax - ymin
  plt.xlim([xmin - 0.05 * xrng, xmax + 0.05 * xrng])
  plt.ylim([ymin - 0.05 * yrng, ymax + 0.05 * yrng])
  plt.savefig(plot_path+'vary_af_expr'+str(expr)+'.pdf')
    
def get_data_vary_af_threshold(input_fn, thresholds, birth, expr):
  """ Vary allele frequency threshold and plot change in identifiability. """
  if birth:
    config.birth_type = birth
    dobs = utils.synthetic_dob(2504)
  else:
    dobs = []
  return get_idabs_vary_af_threshold(input_fn, thresholds, birth, expr, dobs)

def get_idabs_vary_af_threshold(input_fn, thresholds, birth, expr, dobs):
  # Find identifiability of SNPs selected using each threshold.
  birth_type = birth if birth else 'nobirth'
  expr_path = 'vary_af_expr' + str(expr) + '/'
  expr_data_path = data_path + expr_path + 'vary_af_thresh_' + birth_type + '_expr' + str(expr) + '.pickle'
  if os.path.exists(expr_data_path):
    idabs, snp_nums = pickle.load(open(expr_data_path, 'rb'))

  else:
    idabs = []
    snp_nums = []
    for t in thresholds:
      # Extract SNPs according to expriment
      id_fn = ''
      if expr == 1:
        # Extract with AF threshold
        id_fn = af_thresh_fn = input_fn[:-4] + '_af' + str(t) + '.vcf'
        utils.extract_high_af(input_fn, af_thresh_fn, t)
      elif expr == 2:
        # Extract with AF threshold, then extract with SNAP
        af_thresh_fn = input_fn[:-4] + '_af' + str(t) + '.vcf'
        id_fn = redund_fn = af_thresh_fn[:-4] + '_snap.vcf'
        utils.extract_high_af(input_fn, af_thresh_fn, t)
        utils.extract_low_redund(af_thresh_fn, redund_fn, config.snapout_fn)
      elif expr == 3:
        # Extract with SNAP, then extract with AF threshold
        redund_fn = input_fn[:-4] + '_snap.vcf'
        id_fn = af_thresh_fn = redund_fn[:-4] + '_af' + str(t) + '.vcf'
        utils.extract_low_redund(input_fn, redund_fn, config.snapout_fn)
        utils.extract_high_af(redund_fn, af_thresh_fn, t)
      #idab, snp_num = snp_info.identifiability(af_thresh_fn, use_dob=birth, dobs=dobs)

      # Get identifiability data
      idab, snp_num = snp_info.identifiability(id_fn, use_dob=birth, dobs=dobs)
      idabs.append(idab)
      snp_nums.append(snp_num)
    pickle.dump((idabs, snp_nums), open(expr_data_path, 'wb'))

  return idabs, snp_nums

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
  config.snapout_fn = args.snap_output

  if args.snp_sampling:
    snp_sampling(infile)
  if args.af_plot:
    plot_af_per_snp(infile)
  if args.crypto_hash:
    make_hash(infile, outfile)
  if args.af_threshold:
    expr_vary_af_treshold(infile, args.af_threshold)
  if args.entropy:
    expr_top_k_entropy(infile)
  if args.independence:
    expr_vary_rsq(infile)
  if args.entind:
    expr_top_k_independent(infile)

