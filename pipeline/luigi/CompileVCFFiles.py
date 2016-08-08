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



class ConvertLatestClinvarToVCF(luigi.Task):

    def run(self):
      luigi_dir = os.getcwd()

      # Change to correct directory for storing files
      # TODO: maybe we can use temp files or keep in current dir and cleanup after?
      clinvar_file_dir = "../../../brca/pipeline-data/data/ClinVar"
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


# class DownloadLatestESPData(luigi.Task):

#     def run(self):
#       luigi_dir = os.getcwd()

#       esp_file_dir =  "../../../brca/pipeline-data/data/ESP"
#       pipeline_input_dir = "../../../brca/pipeline-data/data/pipeline_input"
#       os.chdir(esp_file_dir)

#       # Download latest gzipped ClinVarFullRelease
#       url = "http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.tar.gz"
#       file_name = url.split('/')[-1]
#       u = urllib2.urlopen(url)
#       f = open(file_name, 'wb')
#       meta = u.info()
#       file_size = int(meta.getheaders("Content-Length")[0])
#       print "Downloading: %s Bytes: %s" % (file_name, file_size)

#       file_size_dl = 0
#       block_sz = 8192
#       while True:
#           buffer = u.read(block_sz)
#           if not buffer:
#               break

#           file_size_dl += len(buffer)
#           f.write(buffer)
#           status = r"%10d  [%3.2f%%]" % (file_size_dl, file_size_dl * 100. / file_size)
#           status = status + chr(8)*(len(status)+1)
#           print status,

#       f.close()
#       print "Finished downloading %s" % (file_name)

#       # Switch back to this file's directory
#       os.chdir(luigi_dir)

class DownloadAndExtractFilesFromESPTar(luigi.Task):

    # NOTE: This task requires some setup to run properly.
    # 1. vcf module must be installed (`pip install PyVCF`)
    # 2. VCFtools must be installed: https://vcftools.github.io/index.html 
    # VCFtools installation is tricky and will require some extra work

    # TODO: break up into smaller pieces and require each as steps
    def run(self):
      luigi_dir = os.getcwd()

      esp_file_dir =  "../../../brca/pipeline-data/data/ESP"
      pipeline_input_dir = "../../../brca/pipeline-data/data/pipeline_input"
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
      os.chdir(luigi_dir)
      esp_method_dir = "../esp"
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


