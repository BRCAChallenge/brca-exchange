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

def uniqueness(input_fn, use_dob, dobs=None):
  """ Computes uniqueness of individuals in the dataset. """
  vcf_reader = vcf.Reader(open(data_path+input_fn, 'rb'))
  m = len(vcf_reader.samples)
  n = sum(1 for _ in vcf_reader)
  vcf_reader = vcf.Reader(open(data_path+input_fn, 'rb'))
  matrix1 = np.matrix(np.zeros(shape=(n, m)), dtype=str)
  matrix2 = np.matrix(np.zeros(shape=(n, m)), dtype=str)

  for i, snp in enumerate(vcf_reader):
    for j, ind in enumerate(snp.samples):
      matrix1[i, j], matrix2[i, j] = ind.gt_bases.split('|')
  matrix1 = matrix1.transpose()
  matrix2 = matrix2.transpose()

  seqs = []
  for i in range(m):
    snp_str1 = ''.join(matrix1[i].tolist()[0])
    snp_str2 = ''.join(matrix2[i].tolist()[0])
    if use_dob:
      seqs.append(snp_str1 + snp_str2 + str(dobs[i]))
    else:
      seqs.append(snp_str1 + snp_str2)
  counter = Counter(seqs)

  o = '' if use_dob else 'o'
  print 'num of unique SNP strings merged, w/' + o + ' DOB', len(counter)
  #print 'num of groups of persons sharing, w/' + o + ' DOB', len([v for v in counter.values() if v > 1])

  return 1.0*len(counter)/m

