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

data_path = config.data_path

def shannon_entropy(reader):
  """ Computes the Shannon Entropy of each SNP.
      Returns a list of entropies.
  """
  entropies = []
  count = 0
  for snp in reader:
    if count > 10:
      #break
      pass
    R = len(snp.alleles)
    probs = utils.allele_frequencies(snp)
    for i in range(R):
      if probs[i] <= 0:
        print snp
    entropy = -1.0 * sum([probs[i] * math.log(probs[i]) if probs[i] > 0 else 0 for i in range(R)])
    entropies.append((count, entropy))
    count += 1
  return entropies

def simpson_index(reader):
  simpson_indices = []
  count = 0
  for snp in reader:
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

def genotypes(input_fn, m, n, use_dob, dobs=None):
  """ Returns concatenated genotype strings
      (and birth information)
      of individuals in input file.
  """
  vcf_reader = vcf.Reader(open(data_path+input_fn, 'rb'))

  # Collect DNA sequences of individuals.
  # There are 2 DNA sequences per person.
  matrix1 = np.matrix(np.zeros(shape=(n, m)), dtype=str)
  matrix2 = np.matrix(np.zeros(shape=(n, m)), dtype=str)
  for i, snp in enumerate(vcf_reader):
    for j, ind in enumerate(snp.samples):
      matrix1[i, j], matrix2[i, j] = ind.gt_bases.split('|')
  matrix1 = matrix1.transpose()
  matrix2 = matrix2.transpose()

  # Genotype sequence of each individual.
  seqs = []
  for i in range(m):
    snp_str1 = ''.join(matrix1[i].tolist()[0])
    snp_str2 = ''.join(matrix2[i].tolist()[0])
    if use_dob:
      seqs.append(snp_str1 + snp_str2 + str(dobs[i]))
    else:
      seqs.append(snp_str1 + snp_str2)
  return seqs

def identifiability(input_fn, use_dob, dobs=None):
  """ Computes identifiability of individuals in the dataset.
      Returns identifiability (proportion of uniquely identified individuals)
      and number of SNPs in input file.
  """
  vcf_reader = vcf.Reader(open(data_path+input_fn, 'rb'))
  m = len(vcf_reader.samples)
  n = sum(1 for _ in vcf_reader)
  seqs = genotypes(input_fn, m, n, use_dob, dobs)
  counter = Counter(seqs)

  o = '' if use_dob else 'o'
  print 'num of unique SNP strings merged, w/' + o + ' DOB', len(counter)
  #print 'num of groups of persons sharing, w/' + o + ' DOB',\
  #      len([v for v in counter.values() if v > 1])

  return 1.0*len(counter)/m, n

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
