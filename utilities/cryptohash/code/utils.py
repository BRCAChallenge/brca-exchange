import vcf
import random
from collections import Counter
import random
import datetime
from datetime import date, timedelta
import calendar
import pickle
import numpy as np
import config
import os
import snap_reader
import time
import snp_info

data_path = config.data_path

def vcf_to_dicts(input_vcf):
  vcf_reader = vcf.Reader(open(data_path + input_vcf, 'rb'))
  snps = []
  for snp in vcf_reader:
    snp_dict = {'snp': snp.ID, samples: []}
  return snps

def extract_low_redund_2(input_fn, snapout_fn, k=-1):
  """ From the input VCF file, return the (top k, if k >= 0) most informative
      rsq-independent SNPs, using given SNAP output file.
  """

  # If k > 0, rank SNPs by entropy.
  if k >= 0:
    ranked_path = input_fn[:-4] + '_ranked.pickle'
    if os.path.exists(data_path + ranked_path):
      snps = pickle.load(open(data_path + ranked_path, 'rb'))
    else:
      vcf_reader = vcf.Reader(open(data_path + input_fn, 'rb'))
      _, snps = snp_info.shannon_entropy(vcf_reader)
      pickle.dump(snps, open(data_path + ranked_path, 'wb'))

  # Else just take SNPs in the order they come.
  else:
    unranked_path = input_fn[:-4] + '_unranked.pickle'
    if os.path.exists(data_path + unranked_path):
      snps = pickle.load(open(data_path + unranked_path, 'rb'))
    else:
      vcf_reader = vcf.Reader(open(data_path + input_fn, 'rb'))
      snps = [snp.ID for snp in vcf_reader]
      pickle.dump(snps, open(data_path + unranked_path, 'wb'))

  # Starting from the highest entropy to the lowest,
  # select each SNP that is not correlated with already selected ones.
  selection = list()
  correlated = set()
  proxy_dict = snap_reader.read_snapout(snapout_fn, input_vcf=input_fn)
  for snp in snps:
    if k > 0 and len(selection) == k:
      break
    if snp in correlated:
      continue
    if snp in proxy_dict:
      c = [proxy['proxy'] for proxy in proxy_dict[snp]]
    else:
      c = []
    correlated.update(c)
    selection.append(snp)

  return selection

def extract_low_redund(input_fn, output_fn, snapout_fn):
  """ Output file with non-redundant SNPs using given SNAP output file. """
  if os.path.exists(data_path + output_fn):
    print output_fn, "exists, not recreating"
    return
  print "creating", output_fn

  proxy_dict = snap_reader.read_snapout(snapout_fn, input_vcf=input_fn)
  vcf_reader = vcf.Reader(open(data_path + input_fn, 'rb'))

  selection = set()
  correlated = set()
  for snp in vcf_reader:
    if snp.ID in correlated:
      continue
    if snp.ID in proxy_dict:
      c = [proxy['proxy'] for proxy in proxy_dict[snp.ID]]
    else:
      c = []
    correlated.update(c)
    selection.add(snp.ID)

  with open(data_path+output_fn, 'wb') as outfile:
    with open(data_path+input_fn, 'rb') as infile:
      for line in infile:
        if line.startswith('#'):
          outfile.write(line)
        else:
          snp_id = line.split('\t')[2]
          if snp_id in selection:
            outfile.write(line)
          else:
            continue

def extract_high_af(input_fn, output_fn, t):
  """ Output file with all SNPs with allele frequencies all >= t% """
  if os.path.exists(data_path + output_fn):
    #print "Regenerating", data_path + output_fn
    print output_fn, "exists, not recreating"
    return

  vcf_reader = vcf.Reader(open(data_path + input_fn, 'rb'))
  snp_indices = []
  for snp_ind, snp in enumerate(vcf_reader):
    if all(af >= 0.01 * t for af in allele_frequencies(snp)):
      snp_indices.append(snp_ind)

  with open(data_path + output_fn, 'wb') as outfile:
    with open(data_path + input_fn, 'rb') as infile:
      snp_ind = 0
      for line in infile:
        if line.startswith('#'):
          outfile.write(line)
        else:
          if snp_ind in snp_indices:
            outfile.write(line)
          snp_ind += 1

def exon_starts_ends(input_fn, names):
  exon_starts = []
  exon_ends = []
  with open(data_path+input_fn, 'rb') as infile:
    for line in infile:
      if line.startswith('#'):
        continue
      fields = line.split('\t')
      name = fields[1]
      if name in names:
        exon_starts += fields[9].split(',')
        exon_ends   += fields[10].split(',')
  assert len(exon_starts) == len(exon_ends), 'exon starts and ends differ in number'
  return exon_starts, exon_ends

exon_starts, exon_ends = exon_starts_ends('RefSeq.genepred', ['NM_000059.3', 'NM_007294.3'])

def extract_exon(input_fn, output_fn):
  """ Extract SNPs on exons and write to new file. """
  if os.path.exists(data_path+output_fn):
    print "Regenerating", data_path+output_fn

  def on_exon(line):
    pos = line.split('\t')[1]
    for i in range(len(exon_starts)):
      start, end = exon_starts[i], exon_ends[i]
      if (pos > start) and (pos < end):
        return True
    return False

  with open(data_path+output_fn, 'wb') as outfile:
    with open(data_path+input_fn, 'rb') as infile:
      for line in infile:
        if (line.startswith('#')) or (on_exon(line)):
          outfile.write(line)

def dump_records(input_fn):
  """ Dump a VCF file as list of vcf.model._Record. """
  vcf_reader = vcf.Reader(open(data_path+input_fn, 'rb'))
  records = [record for record in vcf_reader]
  pickle.dump(records, open(data_path+input_fn+'.pickle', 'wb'))

def make_toy(input_fn, output_fn, p):
  """ Sample lines with p% probability. """
  if os.path.exists(data_path+output_fn):
    print "Regenerating", data_path+output_fn

  snp_num = 0
  snp_sam = 0
  with open(data_path+output_fn, 'wb') as outfile:
    with open(data_path+input_fn, 'rb') as infile:
      for line in infile:
        if line.startswith('#'):
          outfile.write(line)
        else:
          snp_num += 1
          if random.random() < 0.01*p:
            outfile.write(line)
            snp_sam += 1

def allele_frequencies(record):
  """ Allele frequencies of record; [0] is frequency of REF allele. """
  allele_counts = Counter()
  num_chroms = 0.0
  for s in record.samples:
    if s.gt_type is not None:
      for a in s.gt_alleles:
        allele_counts.update([a])
        num_chroms += 1
  return [allele_counts[str(i)]/num_chroms for i in range(0, len(record.ALT)+1)]

def synthetic_headers(line, p):
  """ Generate tab separated headers of synthetic samples. """
  line = line.split('FORMAT')[1]
  sample_num = len(line.split('\t'))
  synth_num = int(sample_num * p * 0.01)
  synth = ['SY'+str(i).zfill(5) for i in range(synth_num)]
  choices = random.sample(xrange(sample_num), synth_num)
  return sample_num, '\t'.join(synth), choices

def synthetic_records(line, person_choices, synthetic):
  """ Generate tab separated SNP values of synthetic coreferent individuals. """
  if person_choices == None:
    return
  fields = line.split('\t')
  alleles = set([fields[3]] + fields[4].split(','))
  sample = line.split('\t')[9:]
  values = [sample[c] for c in person_choices]
  if synthetic: # randomly pick one from possible alleles?
    values = [random.choice(list(alleles - set([values[c]]))) for c in person_choices]
  return '\t'.join(values)

def make_synthetic_coreference(input_fn, output_fn, p):
  """ Generate VCF file with synthetic coreferent records,
      using p% of data in input_fn. """
  if os.path.exists(data_path+output_fn):
    print "Regenerating", data_path+output_fn

  snp_num = 0
  with open(data_path+input_fn, 'rb') as infile:
    for line in infile:
      if not line.startswith('#'):
        snp_num += 1
  snp_choices = random.sample(xrange(snp_num), int(snp_num * p * 0.01))
  print snp_choices
  person_choices = None
  snp_count = 0
  with open(data_path+output_fn, 'wb') as outfile:
    with open(data_path+input_fn, 'rb') as infile:
      for line in infile:
        if line.startswith('##'): # meta information
          outfile.write(line)
        elif line.startswith('#'): # header
          sample_num, syn_headers, person_choices = synthetic_headers(line, 10)
          outfile.write(line.replace('\n', '\t'))
          outfile.write(syn_headers + '\n')
        else:
          outfile.write(line.replace('\n', '\t'))
          outfile.write(synthetic_records(line, person_choices, snp_count in snp_choices) + '\n')
        snp_count += 1
  return sample_num

def all_possible_dobs(min_age, max_age, birth_type):
  """ Return all possible dates of birth of each year, given a range of ages. """
  now = datetime.datetime.now()
  min_year = now.year - max_age
  max_year = now.year - min_age
  births = []
  for year in range(min_year, max_year+1):
    if birth_type == 'date':
      frst_day = date(year, 1, 1)
      last_day = date(year, 12, 31)
      delta = last_day - frst_day
      births += [frst_day + timedelta(days=i) for i in range(delta.days + 1)]
    elif birth_type == 'year':
      births += [year]
  return births

def synthetic_dob(m):
  """ Sample dates of birth of m persons. """
  dobs = all_possible_dobs(20, 65, config.birth_type)
  return np.random.choice(dobs, m, replace=True)


def fill_zeros(m, max_len):
  """ Append zeros to each list in m until they have length max_len;
      return a matrix as an np.array.
  """
  for l in m:
    if len(l) < max_len:
      for i in range(max_len - len(l)):
        l.append(0)
  return np.array([l for l in m])

