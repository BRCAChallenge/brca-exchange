import vcf
import math
from collections import Counter
from pprint import pprint
import numpy as np
import pickle
from matplotlib import pyplot as plt
import utils
import config
import os
from datetime import datetime
import time

data_path = config.data_path

def shannon_entropy(vcf_reader):
  """ Computes the Shannon Entropy of each SNP.
      Returns a list of entropies and a list of RSIDs in corresponding order.
  """
  ents_rsids = []
  #rsids = []
  for snp in vcf_reader:
    R = len(snp.alleles)
    probs = utils.allele_frequencies(snp)
    #for i in range(R):
    #  assert probs[i] >= 0
    entropy = -1.0 * sum([probs[i] * math.log(probs[i]) if probs[i] > 0 else 0 for i in range(R)])
    ents_rsids.append((entropy, snp.ID))
  ents_rsids = sorted(ents_rsids, reverse=True, key=lambda i : i[0])
  entropies = [ent for ent, rsid in ents_rsids]
  rsids = [rsid for ent, rsid in ents_rsids]
  return entropies, rsids

def simpson_index(vcf_reader):
  simpson_indices = []
  count = 0
  for snp in vcf_reader:
    if count > 10:
      #break
      pass
    R = len(snp.alleles)
    probs = utils.allele_frequencies(snp)
    simpson_index = sum([probs[i]**2 for i in range(R)])
    simpson_indices.append((count, simpson_index))
    count += 1
  return simpson_indices

def shannon_and_simpson(use_toy):
  input_fn = '1000G_brca.hg38.vcf'
  toy_fn = input_fn[:-4]+'.toy.vcf' if use_toy else input_fn
  if not os.path.exists(data_path+toy_fn):
    utils.make_toy(input_fn, toy_fn, 20)

  vcf_reader = vcf.Reader(open(data_path+toy_fn, 'rb'))
  shannon = sorted(shannon_entropy(vcf_reader), key=lambda i : i[1], reverse=True)
  vcf_reader = vcf.Reader(open(data_path+toy_fn, 'rb'))
  simpson = sorted(simpson_index(vcf_reader), key=lambda i : i[1], reverse=False)
  pprint(zip(shannon, simpson))
  l1 = [shannon[i][0] for i in range(len(shannon))]
  l2 = [simpson[i][0] for i in range(len(simpson))]
  #print l1
  #print l2
  #print zip(l1, l2)
  #pprint(shannon)
  #pprint(simpson)
  # => Entropy and Simpson give similar ranking; expected as both based on probs

def genotypes(input_fn, m, n, use_dob, dobs=None, snps=None):
  """ Returns concatenated genotype strings
      (and birth information)
      of individuals in input file.
  """
  genotype_path = input_fn[:-4] + '_genotype.pickle'
  # If not given set of SNPs to generate genotypes, use saved genotypes.
  if snps==None and os.path.exists(genotype_path):
    return pickle.load(open(data_path + genotype_path, 'rb'))

  vcf_reader = vcf.Reader(open(data_path+input_fn, 'rb'))
  # Collect DNA sequences of individuals.
  # There are 2 DNA sequences per person.
  matrix1 = np.matrix(np.zeros(shape=(m, n)), dtype=str)
  matrix2 = np.matrix(np.zeros(shape=(m, n)), dtype=str)
  for i, snp in enumerate(vcf_reader):
    if snps==None or snp.ID in snps: # snps=None or snp.ID in snps
      for j, ind in enumerate(snp.samples):
        matrix1[j, i], matrix2[j, i] = ind.gt_bases.split('|')

  # Genotype sequence of each individual.
  seqs = []
  for i in range(m):
    snp_str1 = ''.join(matrix1[i].tolist()[0])
    snp_str2 = ''.join(matrix2[i].tolist()[0])
    if dobs != None and len(dobs)>0:
      seqs.append(snp_str1 + snp_str2 + str(dobs[i]))
    else:
      seqs.append(snp_str1 + snp_str2)
  pickle.dump(seqs, open(data_path + genotype_path, 'wb'))
  return seqs

def identifiability(input_fn, use_dob, dobs=None, snps=None):
  """ Computes identifiability of individuals in the dataset.
      Returns identifiability (proportion of uniquely identified individuals)
      and number of SNPs in input file.
  """
  start_time = time.time()
  mn_path = input_fn[:-4] + '_mn.pickle'
  if os.path.exists(data_path + mn_path):
    m, n = pickle.load(open(data_path + mn_path, 'rb'))
  else:
    vcf_reader = vcf.Reader(open(data_path + input_fn, 'rb'))
    m = len(vcf_reader.samples)
    n = sum(1 for _ in vcf_reader)
    pickle.dump((m, n), open(data_path + mn_path, 'wb'))
  seqs = genotypes(input_fn, m, n, use_dob, dobs, snps)
  counter = Counter(seqs)

  o = '' if use_dob else 'o'
  print 'num of unique SNP strings merged, w/' + o + ' DOB', len(counter)
  #print 'num of groups of persons sharing, w/' + o + ' DOB',\
  #      len([v for v in counter.values() if v > 1])

  print 'identifiability', time.time() - start_time
  return 1.0*len(counter)/m, len(snps) if snps else n

def score_snp(snp):
  """ Score a SNP on its identifiability of individuals.
      Returns a double score s, 0 <= s <= 1.
  """
  pass

def conditional(input_fn):
  """ Computes conditional distribution of each SNP:
      how informative the allele on one chromotid is
      about the allele on the other chromotid.
  """
  vcf_reader = vcf.Reader(open(data_path+input_fn, 'rb'))
  dists_a = []
  dists_b = []
  dists_ab = []
  num_snp = sum([1 for _ in vcf_reader])
  vcf_reader = vcf.Reader(open(data_path+input_fn, 'rb'))
  
  #                            ?
  # For each SNP, for each chromotid, probability of each observed value,
  #           and for both chromotids, joint probability
  prior_dists_file = 'prior_dists_a_b_ab.pickle'
  if not os.path.exists(data_path + prior_dists_file):
    for snp in vcf_reader:
      hist_a = {}
      hist_b = {}
      hist_ab = {}
      for ind in snp.samples:
        bases = ind.gt_bases.split('|')
        a = bases[0]
        b = bases[1]
        ab = tuple(bases)
        hist_a.update({a : hist_a.get(a, 0) + 1})
        hist_b.update({b : hist_b.get(b, 0) + 1})
        hist_ab.update({ab : hist_ab.get(ab, 0) + 1})
      dists_a.append(hist_a)
      dists_b.append(hist_b)
      dists_ab.append(hist_ab)
    prior_dists = (dists_a, dists_b, dists_ab)
    pickle.dump(prior_dists, open(data_path + prior_dists_file, 'wb'))
  
  else:
    prior_dists = pickle.load(open(data_path + prior_dists_file, 'rb'))

  for dists in prior_dists:
    for hist in dists:
      v_sum = sum(hist.values())
      for k, v in hist.items():
        hist.update({k : 1.0*v/v_sum})
  print 'done'

  dists_a, dists_b, dists_ab = prior_dists
  
  def get_conditionals(joints, margs, num_snp, marg_ind):
    assert marg_ind == 0 or marg_ind == 1, 'Bad index of marginal'
    if marg_ind == 0:
      cond_file = 'cond_dists_b_a.pickle'
    else:
      cond_file = 'cond_dists_a_b.pickle'

    if not os.path.exists(data_path + cond_file):
      conditionals = []
      for i in range(num_snp):
        joint = joints[i]
        marg = margs[i]
        cond = {}
        for m in marg.keys():
          for j in joint.keys():
            if j[marg_ind] == m:
              if marg_ind == 0:
                c = j[1]+'_given_'+m
              else:
                c = m+'_given_'+j[0]
              cond.update({c : joint[j] / marg[m]})
        print cond
        conditionals.append(cond)
        pickle.dump(conditionals, open(data_path + cond_file, 'wb'))
    else:
      conditionals = pickle.load(open(data_path + cond_file, 'rb'))
    return conditionals

  # P(b|a) = P(ab) / P(a)
  b_given_as = get_conditionals(dists_ab, dists_a, num_snp, 0)
  # P(a|b) = P(ab) / P(b)
  a_given_bs = get_conditionals(dists_ab, dists_b, num_snp, 1)

def ld_information(input_fn, proxy_dict):
  """ Compute the LD score. """
  def score_ld(rsid, proxy_dict):
    proxies = proxy_dict[rsid]
    num_prx = len(proxies)
    sum_rsq = 0
    for proxy in proxies:
      rsq = float(proxy['rsquared'])
      sum_rsq += rsq
    mean_rsq = sum_rsq/num_prx
    return mean_rsq

  vcf_reader = vcf.Reader(open(data_path+input_fn, 'rb'))
  found = 0
  snp_ld_dict = {}
  for snp in vcf_reader:
    if snp.ID in proxy_dict.keys():
      found += 1
      score = score_ld(snp.ID, proxy_dict)
      #print score
    else:
      score = 0
    snp_ld_dict[snp.ID] = score
  #pprint(snp_ld_dict)
  found = 1.0 * sum([v for v in snp_ld_dict.values() if v != 0])
  print 100*found/len(snp_ld_dict), '% of SNPs in dataset found in SNAP output'
  print [k for k, v in snp_ld_dict.items() if v != 0]
  return [k for k, v in snp_ld_dict.items() if v != 0]
