import subprocess
import os
import urllib2
import tarfile
import luigi

from shutil import copyfile

# Import methods from pipeline files
import sys
sys.path.insert(0, '../esp')
import espExtract

# Environment variables
# TODO: Formalize for other systems
# NOTE: These must be changed to represent the directory structure of the machine running this script
os.environ['BRCA_RESOURCES'] = '/Users/zackfischmann/UCSC/brca/brca-resources'
os.environ['PIPELINE_INPUT'] = '/Users/zackfischmann/UCSC/brca/pipeline-data/data/pipeline_input'
os.environ['BIC'] = '/Users/zackfischmann/UCSC/brca/pipeline-data/data/BIC'
os.environ['CLINVAR'] = '/Users/zackfischmann/UCSC/brca/pipeline-data/data/ClinVar'
os.environ['ESP'] = '/Users/zackfischmann/UCSC/brca/pipeline-data/data/ESP'
os.environ['LUIGI'] = '/Users/zackfischmann/UCSC/brca-exchange/pipeline/luigi'
os.environ['BIC_METHODS'] = '/Users/zackfischmann/UCSC/brca-exchange/pipeline/bic'
os.environ['ESP_METHODS'] = '/Users/zackfischmann/UCSC/brca-exchange/pipeline/esp'

# Globals
brca_resources_dir = os.environ['BRCA_RESOURCES']
pipeline_input_dir = os.environ['PIPELINE_INPUT']
bic_file_dir =  os.environ['BIC']
clinvar_file_dir = os.environ['CLINVAR']
esp_file_dir = os.environ['ESP']
luigi_dir = os.environ['LUIGI']
bic_method_dir = os.environ['BIC_METHODS']
esp_method_dir = os.environ['ESP_METHODS']

class ConvertLatestClinvarToVCF(luigi.Task):

    def run(self):
      # Change to correct directory for storing files
      # TODO: maybe we can use temp files or keep in current dir and cleanup after?
      os.chdir(clinvar_file_dir)

      # Download latest gzipped ClinVarFullRelease
      url = "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz"
      file_name = url.split('/')[-1]
      u = urllib2.urlopen(url)
      f = open(file_name, 'wb')
      meta = u.info()
      file_size = int(meta.getheaders("Content-Length")[0])
      print "Downloading: %s Bytes: %s" % (file_name, file_size)

      file_size_dl = 0
      block_sz = 8192
      while True:
          buffer = u.read(block_sz)
          if not buffer:
              break

          file_size_dl += len(buffer)
          f.write(buffer)
          status = r"%10d  [%3.2f%%]" % (file_size_dl, file_size_dl * 100. / file_size)
          status = status + chr(8)*(len(status)+1)
          print status,

      f.close()
      print "Finished downloading %s" % (file_name)

      # Switch back to this file's directory
      os.chdir(luigi_dir)

      # Convert gzipped file to vcf using makefile in pipeline/clinvar/
      print "Converting %s to VCF format... this takes a while." % (file_name)
      sp = subprocess.Popen(["make"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd="../clinvar")
      out, err = sp.communicate()
      if out:
          print "standard output of subprocess:"
          print out
      if err:
          print "standard error of subprocess:"
          print err


class DownloadAndExtractFilesFromESPTar(luigi.Task):

    # NOTE: This task requires some setup to run properly.
    # 1. vcf module must be installed (`pip install PyVCF`)
    # 2. VCFtools must be installed: https://vcftools.github.io/index.html 
    # VCFtools installation is tricky and will require some extra work

    # TODO: break up into smaller pieces and require each as steps
    def run(self):
      os.chdir(esp_file_dir)

      url = "http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.tar.gz"
      file_name = url.split('/')[-1]

      u = urllib2.urlopen(url)
      f = open(file_name, 'wb')
      meta = u.info()
      file_size = int(meta.getheaders("Content-Length")[0])
      print "Downloading: %s Bytes: %s" % (file_name, file_size)

      file_size_dl = 0
      block_sz = 8192
      while True:
          buffer = u.read(block_sz)
          if not buffer:
              break

          file_size_dl += len(buffer)
          f.write(buffer)
          status = r"%10d  [%3.2f%%]" % (file_size_dl, file_size_dl * 100. / file_size)
          status = status + chr(8)*(len(status)+1)
          print status,

      f.close()
      print "Finished downloading %s" % (file_name)

      # Extract contents of tarfile
      tar = tarfile.open(file_name, "r:gz")
      tar.extractall()
      tar.close()
      print "Finished extracting files from %s" % (file_name)

      # Switch back to this file's directory
      os.chdir(esp_method_dir)

      # Extract data for BRCA1 region
      brca1_region_file = esp_file_dir + "/ESP6500SI-V2-SSA137.GRCh38-liftover.chr17.snps_indels.vcf"
      brca1_region_output = esp_file_dir + "/esp.brca1.vcf"
      args = ["python", "espExtract.py", brca1_region_file, "--start", "43044295", "--end", "43125483", "--full", "1", "-o", brca1_region_output]
      print "Calling espExtract.py for BRCA1 region with the following arguments: %s" % (args)
      sp = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = sp.communicate()
      if out:
          print "standard output of subprocess for brca1:"
          print out
      if err:
          print "standard error of subprocess for brca1:"
          print err

      # Extract data for BRCA2 region
      brca2_region_file = esp_file_dir + '/ESP6500SI-V2-SSA137.GRCh38-liftover.chr13.snps_indels.vcf'
      brca2_region_output = esp_file_dir + "/esp.brca2.vcf"
      args = ["python", "espExtract.py", brca2_region_file, "--start", "43044295", "--end", "43125483", "--full", "1", "-o", brca2_region_output]
      print "Calling espExtract.py for BRCA 2 region with the following arguments: %s" % (args)
      sp = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = sp.communicate()
      if out:
          print "standard output of subprocess for brca2:"
          print out
      if err:
          print "standard error of subprocess for brca2:"
          print err

      # Concatenate extracted BRCA1/BRCA2 region data
      # Note: requires correct installation of VCF tools and export PERL5LIB=/path/to/your/vcftools-directory/src/perl/ in path
      concatenated_brca_output_file = esp_file_dir + "/esp.brca12.hg38.vcf"
      writable_concatenated_brca_output_file = open(concatenated_brca_output_file, 'w')
      args = ["vcf-concat", brca1_region_output, brca2_region_output]
      print "Calling vcf-concat with the following args: %s" % (args)
      sp = subprocess.Popen(args, stdout=writable_concatenated_brca_output_file, stderr=subprocess.PIPE)
      out, err = sp.communicate()
      if out:
          print "standard output of subprocess:"
          print out
      if err:
          print "standard error of subprocess:"
          print err

      writable_concatenated_brca_output_file.close()
      print "Concatenation complete."

      # Sort concatenated BRCA1/2 data
      sorted_concatenated_brca_output_file = esp_file_dir + "/esp.brca12.sorted.hg38.vcf"
      writable_sorted_concatenated_brca_output_file = open(sorted_concatenated_brca_output_file, 'w')
      args = ["vcf-sort", concatenated_brca_output_file]
      print "Calling vcf-sort with the following args: %s" % (args)
      sp = subprocess.Popen(args, stdout=writable_sorted_concatenated_brca_output_file, stderr=subprocess.PIPE)
      out, err = sp.communicate()
      if out:
          print "standard output of subprocess:"
          print out
      if err:
          print "standard error of subprocess:"
          print err

      writable_sorted_concatenated_brca_output_file.close()
      print "Sorting of concatenated files complete."

      # Copy sorted data to correct path
      print "Copying sorted and concatenated VCF to the final pipeline_input destination"
      final_destination = pipeline_input_dir + "/esp.brca12.sorted.hg38.vcf"
      copyfile(sorted_concatenated_brca_output_file, final_destination)
      print "File copied to %s" % (final_destination)


class DownloadAndExtractFilesFromBIC(luigi.Task):

    # TODO: Figure out how to store u/p for safe retrieval
    # NOTE: U/P can be found in /hive/groups/cgl/brca/phase1/data/bic/account.txt at UCSC
    username = ''
    password = ''

    def run(self):
      os.chdir(bic_file_dir)

      brca1_data_url = "https://research.nhgri.nih.gov/projects/bic/Member/cgi-bin/bic_query_result.cgi/brca1_data.txt?table=brca1_exons&download=1&submit=Download"
      brca1_file_name = "brca1_data.txt"

      brca2_data_url = "https://research.nhgri.nih.gov/projects/bic/Member/cgi-bin/bic_query_result.cgi/brca2_data.txt?table=brca2_exons&download=1&submit=Download"
      brca2_file_name = "brca2_data.txt"

      # Download brca1 data
      p = urllib2.HTTPPasswordMgrWithDefaultRealm()

      p.add_password(None, brca1_data_url, self.username, self.password)

      handler = urllib2.HTTPBasicAuthHandler(p)
      opener = urllib2.build_opener(handler)
      urllib2.install_opener(opener)

      data = urllib2.urlopen(brca1_data_url).read()
      f = open(brca1_file_name, "wb")
      f.write(data)
      f.close()
      print "Finished downloading %s" % (brca1_file_name)

      # Download brca2 data
      p = urllib2.HTTPPasswordMgrWithDefaultRealm()

      p.add_password(None, brca2_data_url, self.username, self.password)

      handler = urllib2.HTTPBasicAuthHandler(p)
      opener = urllib2.build_opener(handler)
      urllib2.install_opener(opener)

      data = urllib2.urlopen(brca2_data_url).read()
      f = open(brca2_file_name, "wb")
      f.write(data)
      f.close()
      print "Finished downloading %s" % (brca2_file_name)

      # Convert files to vcf
      os.chdir(bic_method_dir)

      bic_vcf_file = bic_file_dir + "/bicSnp.hg19.vcf"
      writable_bic_vcf_file = open(bic_vcf_file, 'w')
      args = ["python", "convBic.py", "--brca1", bic_file_dir + "/" + brca1_file_name, "--brca2", bic_file_dir + "/" + brca2_file_name]
      print "Converting BRCA1/2 BIC data to vcf with the following args: %s" % (args)
      sp = subprocess.Popen(args, stdout=writable_bic_vcf_file, stderr=subprocess.PIPE)
      out, err = sp.communicate()
      if out:
          print "standard output of subprocess:"
          print out
      if err:
          print "standard error of subprocess:"
          print err

      writable_bic_vcf_file.close()
      print "Conversion of BIC data to VCF complete."

      # Note: CrossMap.py must be installed `pip install CrossMap` for the next step
      os.chdir(luigi_dir)
      bic_hg38_vcf_file = bic_file_dir + "/bicSnp.hg38.vcf"
      writable_bic_hg38_vcf_file = open(bic_hg38_vcf_file, 'w')
      args = ["CrossMap.py", "vcf", brca_resources_dir + "/hg19ToHg38.over.chain.gz", bic_vcf_file, brca_resources_dir + "/hg38.fa", bic_hg38_vcf_file]
      print "Running crossmap.py with the following args: %s" % (args)
      sp = subprocess.Popen(args, stdout=writable_bic_hg38_vcf_file, stderr=subprocess.PIPE)
      out, err = sp.communicate()
      if out:
          print "standard output of subprocess:"
          print out
      if err:
          print "standard error of subprocess:"
          print err
      print "Crossmap.py execution complete, created bicSnp.hg38.vcf."

      # Sort hg38 vcf file, requires VCFtools.
      sorted_hg38_vcf_file = bic_file_dir + "/bicSnp.sorted.hg38.vcf"
      writable_sorted_hg38_vcf_file = open(sorted_hg38_vcf_file, 'w')
      args = ["vcf-sort", bic_hg38_vcf_file]
      print "Running vcf-sort with the following args: %s" % (args)
      sp = subprocess.Popen(args, stdout=writable_sorted_hg38_vcf_file, stderr=subprocess.PIPE)
      out, err = sp.communicate()
      if out:
          print "standard output of subprocess:"
          print out
      if err:
          print "standard error of subprocess:"
          print err
      print "Sorting of hg38 vcf file complete."

      # Copy sorted file to pipeline_input directory
      print "Copying sorted hg38 VCF to the final pipeline_input destination"
      final_destination = pipeline_input_dir + "/bicSnp.sorted.hg38.vcf"
      copyfile(sorted_hg38_vcf_file, final_destination)
      print "File copied to %s" % (final_destination)
