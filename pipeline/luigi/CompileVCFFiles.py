import subprocess
import os
import urllib2
import tarfile
import datetime
import luigi

from shutil import copyfile

# Import methods from pipeline files
import sys
sys.path.insert(0, '../esp')
import espExtract

def create_path_if_nonexistent(path):
  if not os.path.exists(path):
    os.makedirs(path)
  return path

# Create directories, environment variables, and globals for scripts to run.
brca_resources_dir = os.environ['BRCA_RESOURCES'] = create_path_if_nonexistent(os.path.abspath('../brca/brca-resources'))
pipeline_input_dir = os.environ['PIPELINE_INPUT'] = create_path_if_nonexistent(os.path.abspath('../brca/pipeline-data/data/pipeline_input'))
brca_pipeline_data_dir = os.environ['BRCA_PIPELINE_DATA'] = create_path_if_nonexistent(os.path.abspath('../brca/pipeline-data'))

bic_file_dir = os.environ['BIC'] = create_path_if_nonexistent(os.path.abspath('../brca/pipeline-data/data/BIC'))
clinvar_file_dir = os.environ['CLINVAR'] = create_path_if_nonexistent(os.path.abspath('../brca/pipeline-data/data/ClinVar'))
esp_file_dir = os.environ['ESP'] = create_path_if_nonexistent(os.path.abspath('../brca/pipeline-data/data/ESP'))
lovd_file_dir = os.environ['LOVD'] = create_path_if_nonexistent(os.path.abspath('../brca/pipeline-data/data/LOVD'))
ex_lovd_file_dir = os.environ['EXLOVD'] = create_path_if_nonexistent(os.path.abspath('../brca/pipeline-data/data/exLOVD'))
exac_file_dir = os.environ['EXAC'] = create_path_if_nonexistent(os.path.abspath('../brca/pipeline-data/data/exac'))
g1k_file_dir = os.environ['G1K'] = create_path_if_nonexistent(os.path.abspath('../brca/pipeline-data/data/G1K'))

luigi_dir = os.environ['LUIGI'] = os.getcwd()

bic_method_dir = os.environ['BIC_METHODS'] = os.path.abspath('../bic')
esp_method_dir = os.environ['ESP_METHODS'] = os.path.abspath('../esp')
lovd_method_dir = os.environ['LOVD_METHODS'] = os.path.abspath('../lovd')
g1k_method_dir = os.environ['G1K_METHODS'] = os.path.abspath('../1000_Genomes')


class ConvertLatestClinvarToVCF(luigi.Task):

    def run(self):
      os.chdir(clinvar_file_dir)
      
      # Create required xml file for make to run properly
      clinvar_xml_file = clinvar_file_dir + "/ClinVarBrca.xml"
      if not os.path.exists(clinvar_xml_file):
        open(clinvar_xml_file, 'w').close() 

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

      # Makefile requires a particular environment variable for output
      my_env = os.environ.copy()
      my_env["BRCA_PIPELINE_DATA"] = brca_pipeline_data_dir
      
      # Convert gzipped file to vcf using makefile in pipeline/clinvar/
      print "Converting %s to VCF format. This takes a while..." % (file_name)
      sp = subprocess.Popen(["make"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd="../clinvar", env=my_env)
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

      # Download ESP data
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
    date = luigi.DateParameter(default=datetime.date.today())
    u = luigi.Parameter()
    p = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(pipeline_input_dir + "/bicSnp.sorted.hg38.vcf")

    def run(self):
      os.chdir(bic_file_dir)

      brca1_data_url = "https://research.nhgri.nih.gov/projects/bic/Member/cgi-bin/bic_query_result.cgi/brca1_data.txt?table=brca1_exons&download=1&submit=Download"
      brca1_file_name = "brca1_data.txt"

      brca2_data_url = "https://research.nhgri.nih.gov/projects/bic/Member/cgi-bin/bic_query_result.cgi/brca2_data.txt?table=brca2_exons&download=1&submit=Download"
      brca2_file_name = "brca2_data.txt"

      # Download brca1 data
      p = urllib2.HTTPPasswordMgrWithDefaultRealm()

      p.add_password(None, brca1_data_url, self.u, self.p)

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

      p.add_password(None, brca2_data_url, self.u, self.p)

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
      with self.output().open("w") as vcf_file:
        args = ["vcf-sort", bic_hg38_vcf_file]
        print "Running vcf-sort with the following args: %s" % (args)
        sp = subprocess.Popen(args, stdout=vcf_file, stderr=subprocess.PIPE)
        out, err = sp.communicate()
        if out:
            print "standard output of subprocess:"
            print out
        if err:
            print "standard error of subprocess:"
            print err
      print "Sorting of hg38 vcf file complete."

class ExtractAndConvertFilesFromEXLOVD(luigi.Task):
    
    def run(self):
      os.chdir(lovd_method_dir)
      
      # extract_data.py -u http://hci-exlovd.hci.utah.edu/ -l BRCA1 BRCA2 -o $EXLOVD
      ex_lovd_data_host_url = "http://hci-exlovd.hci.utah.edu/"
      args = ["extract_data.py", "-u", ex_lovd_data_host_url, "-l", "BRCA1", "BRCA2", "-o", ex_lovd_file_dir]
      print "Running extract_data.py with the following args: %s" % (args)
      sp = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = sp.communicate()
      if out:
          print "standard output of subprocess:"
          print out
      if err:
          print "standard error of subprocess:"
          print err
      print "Extracted data from %s." % (ex_lovd_data_host_url)

      # Convert extracted flat file to vcf format ./lovd2vcf -i $EXLOVD/BRCA1.txt -o $EXLOVD/exLOVD_brca1.hg19.vcf -a exLOVDAnnotation -b 1 -r $BRCA_RESOURCES/refseq_annotation.hg19.gp -g $BRCA_RESOURCES/hg19.fa
      args = ["./lovd2vcf", "-i", ex_lovd_file_dir + "/BRCA1.txt", "-o", ex_lovd_file_dir + "/exLOVD_brca1.hg19.vcf", "-a", "exLOVDAnnotation", "-b", "1", "-r", brca_resources_dir + "/refseq_annotation.hg19.gp", "-g", brca_resources_dir + "/hg19.fa"]
      print "Running lovd2vcf with the following args: %s" % (args)
      sp = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = sp.communicate()
      if out:
          print "standard output of subprocess:"
          print out
      if err:
          print "standard error of subprocess:"
          print err
      print "Converted extracted BRCA1 flat file to vcf format."

      # ./lovd2vcf -i output_directory/BRCA2.txt -o exLOVD_brca2.vcf -a $EXLOVD/exLOVDAnnotation -b 2 -r $BRCA_RESOURCES/refseq_annotation.hg19.gp -g $BRCA_RESOURCES/hg19.fa
      args = ["./lovd2vcf", "-i", ex_lovd_file_dir + "/BRCA2.txt", "-o", "exLOVD_brca2.vcf", "-a", ex_lovd_file_dir + "exLOVDAnnotation", "-b", "2", "-r", brca_resources_dir + "/refseq_annotation.hg19.gp", "-g", brca_resources_dir + "/hg19.fa"]
      print "Running lovd2vcf with the following args: %s" % (args)
      sp = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = sp.communicate()
      if out:
          print "standard output of subprocess:"
          print out
      if err:
          print "standard error of subprocess:"
          print err
      print "Converted extracted BRCA2 flat file to vcf format."

      # vcf-concat $EXLOVD/exLOVD_brca1.hg19.vcf $EXLOVD/exLOVD_brca2.hg19.vcf > $EXLOVD/exLOVD_brca12.hg19.vcf
      ex_lovd_brca12_hg19_vcf_file = ex_lovd_file_dir + "/exLOVD_brca12.hg19.vcf"
      writable_ex_lovd_brca12_hg19_vcf_file = open(ex_lovd_brca12_hg19_vcf_file, 'w')
      args = ["vcf-concat", ex_lovd_file_dir + "/exLOVD_brca1.hg19.vcf", ex_lovd_file_dir + "/exLOVD_brca2.hg19.vcf"]
      print "Running lovd2vcf with the following args: %s" % (args)
      sp = subprocess.Popen(args, stdout=writable_ex_lovd_brca12_hg19_vcf_file, stderr=subprocess.PIPE)
      out, err = sp.communicate()
      if out:
          print "standard output of subprocess:"
          print out
      if err:
          print "standard error of subprocess:"
          print err
      print "Concatenated BRCA1 and BRCA2 vcf files into %s." % (ex_lovd_brca12_hg19_vcf_file)
      
      # `CrossMap.py vcf $BRCA_RESOURCES/hg19ToHg38.over.chain.gz $EXLOVD/exLOVD_brca12.hg19.vcf $BRCA_RESOURCES/hg38.fa $EXLOVD/exLOVD_brca12.hg38.vcf
      args = ["CrossMap.py", "vcf", brca_resources_dir + "/hg19ToHg38.over.chain.gz", ex_lovd_brca12_hg19_vcf_file, brca_resources_dir + "/hg38.fa", ex_lovd_file_dir + "/exLOVD_brca12.hg38.vcf"]
      print "Running CrossMap.py with the following args: %s" % (args)
      sp = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = sp.communicate()
      if out:
          print "standard output of subprocess:"
          print out
      if err:
          print "standard error of subprocess:"
          print err
      print "Crossmap complete."
      
      # vcf-sort $EXLOVD/exLOVD_brca12.hg38.vcf > $EXLOVD/exLOVD_brca12.sorted.hg38.vcf
      sorted_ex_lovd_brca12_hg38_vcf_file = ex_lovd_file_dir + "/exLOVD_brca12.sorted.hg38.vcf"
      writable_sorted_ex_lovd_brca12_hg38_vcf_file = open(sorted_ex_lovd_brca12_hg38_vcf_file, 'w')
      args = ["vcf-sort", ex_lovd_file_dir + "/exLOVD_brca12.hg38.vcf"]
      print "Running lovd2vcf with the following args: %s" % (args)
      sp = subprocess.Popen(args, stdout=writable_sorted_ex_lovd_brca12_hg38_vcf_file, stderr=subprocess.PIPE)
      out, err = sp.communicate()
      if out:
          print "standard output of subprocess:"
          print out
      if err:
          print "standard error of subprocess:"
          print err
      print "Sorted BRCA1/2 hg38 vcf file into %s." % (sorted_ex_lovd_brca12_hg38_vcf_file)
      
      # `cp $EXLOVD/exLOVD_brca12.sorted.hg38.vcf $PIPELINE_INPUT
      print "Copying sorted hg38 VCF to the final pipeline_input destination"
      final_destination = pipeline_input_dir + "/exLOVD_brca12.sorted.hg38.vcf"
      copyfile(sorted_ex_lovd_brca12_hg38_vcf_file, final_destination)
      print "File copied to %s" % (final_destination)

class ExtractAndConvertFilesFromLOVD(luigi.Task):
    
    def run(self):
      os.chdir(lovd_method_dir)
      
      # extract_data.py -u http://databases.lovd.nl/shared/ -l BRCA1 BRCA2 -o $LOVD
      lovd_data_host_url = "http://databases.lovd.nl/shared/"
      args = ["extract_data.py", "-u", lovd_data_host_url, "-l", "BRCA1", "BRCA2", "-o", lovd_file_dir]
      print "Running extract_data.py with the following args: %s" % (args)
      sp = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = sp.communicate()
      if out:
          print "standard output of subprocess:"
          print out
      if err:
          print "standard error of subprocess:"
          print err
      print "Extracted data from %s." % (lovd_data_host_url)

      # ./lovd2vcf -i $LOVD/BRCA1.txt -o $LOVD/sharedLOVD_brca1.hg19.vcf -a sharedLOVDAnnotation -b 1 -r $BRCA_RESOURCES/refseq_annotation.hg19.gp -g $BRCA_RESOURCES/hg19.fa
      args = ["./lovd2vcf", "-i", lovd_file_dir + "/BRCA1.txt", "-o", lovd_file_dir + "/sharedLOVD_brca1.hg19.vcf", "-a", "sharedLOVDAnnotation", "-b", "1", "-r", brca_resources_dir + "/refseq_annotation.hg19.gp", "-g", brca_resources_dir + "/hg19.fa"]
      print "Running lovd2vcf with the following args: %s" % (args)
      sp = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = sp.communicate()
      if out:
          print "standard output of subprocess:"
          print out
      if err:
          print "standard error of subprocess:"
          print err
      print "Converted extracted BRCA1 flat file to vcf format."

      # ./lovd2vcf -i $LOVD/BRCA2.txt -o $LOVD/sharedLOVD_brca2.hg19.vcf -a sharedLOVDAnnotation -b 2 -r $BRCA_RESOURCES/refseq_annotation.hg19.gp -g $BRCA_RESOURCES/hg19.fa 
      # Comment: BRCA2.txt might or might not be created. It seems like there's an error condition that kills BRCA2. It's not worth fixing the error condition at this point.
      args = ["./lovd2vcf", "-i", lovd_file_dir + "/BRCA2.txt", "-o", "sharedLOVD_brca2.hg19.vcf", "-a", "sharedLOVDAnnotation", "-b", "2", "-r", brca_resources_dir + "/refseq_annotation.hg19.gp", "-g", brca_resources_dir + "/hg19.fa"]
      print "Running lovd2vcf with the following args: %s" % (args)
      sp = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = sp.communicate()
      if out:
          print "standard output of subprocess:"
          print out
      if err:
          print "standard error of subprocess:"
          print err
      print "Converted extracted BRCA2 flat file to vcf format."

      # vcf-concat $LOVD/sharedLOVD_brca1.hg19.vcf $LOVD/sharedLOVD_brca2.hg19.vcf > $LOVD/sharedLOVD_brca12.hg19.vcf
      shared_lovd_brca12_hg19_vcf_file = lovd_file_dir + "/sharedLOVD_brca12.hg19.vcf"
      writable_shared_lovd_brca12_hg19_vcf_file = open(shared_lovd_brca12_hg19_vcf_file, 'w')
      args = ["vcf-concat", lovd_file_dir + "/sharedLOVD_brca1.hg19.vcf", lovd_file_dir + "/sharedLOVD_brca2.hg19.vcf"]
      print "Running lovd2vcf with the following args: %s" % (args)
      sp = subprocess.Popen(args, stdout=writable_shared_lovd_brca12_hg19_vcf_file, stderr=subprocess.PIPE)
      out, err = sp.communicate()
      if out:
          print "standard output of subprocess:"
          print out
      if err:
          print "standard error of subprocess:"
          print err
      print "Concatenated BRCA1 and BRCA2 vcf files into %s." % (shared_lovd_brca12_hg19_vcf_file)
      
      # `CrossMap.py vcf $BRCA_RESOURCES/hg19ToHg38.over.chain.gz $LOVD/sharedLOVD_brca12.hg19.vcf $BRCA_RESOURCES/hg38.fa $LOVD/sharedLOVD_brca12.hg38.vcf
      args = ["CrossMap.py", "vcf", brca_resources_dir + "/hg19ToHg38.over.chain.gz", shared_lovd_brca12_hg19_vcf_file, brca_resources_dir + "/hg38.fa", lovd_file_dir + "/sharedLOVD_brca12.hg38.vcf"]
      print "Running CrossMap.py with the following args: %s" % (args)
      sp = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = sp.communicate()
      if out:
          print "standard output of subprocess:"
          print out
      if err:
          print "standard error of subprocess:"
          print err
      print "Crossmap complete."
      
      # vcf-sort $LOVD/sharedLOVD_brca12.38.vcf > $LOVD/sharedLOVD_brca12.sorted.38.vcf
      sorted_shared_lovd_brca12_hg38_vcf_file = lovd_file_dir + "/sharedLOVD_brca12.sorted.hg38.vcf"
      writable_sorted_shared_lovd_brca12_hg38_vcf_file = open(sorted_shared_lovd_brca12_hg38_vcf_file, 'w')
      args = ["vcf-sort", lovd_file_dir + "/sharedLOVD_brca12.hg38.vcf"]
      print "Running lovd2vcf with the following args: %s" % (args)
      sp = subprocess.Popen(args, stdout=writable_sorted_shared_lovd_brca12_hg38_vcf_file, stderr=subprocess.PIPE)
      out, err = sp.communicate()
      if out:
          print "standard output of subprocess:"
          print out
      if err:
          print "standard error of subprocess:"
          print err
      print "Sorted BRCA1/2 hg38 vcf file into %s." % (writable_sorted_shared_lovd_brca12_hg38_vcf_file)
      
      # `cp $LOVD/sharedLOVD_brca12.sorted.38.vcf $PIPELINE_INPUT
      print "Copying sorted hg38 VCF to the final pipeline_input destination"
      final_destination = pipeline_input_dir + "/sharedLOVD_brca12.sorted.hg38.vcf"
      copyfile(sorted_shared_lovd_brca12_hg38_vcf_file, final_destination)
      print "File copied to %s" % (final_destination)

class DownloadAndExtractFilesFromG1K(luigi.Task):

    def run(self):
      os.chdir(g1k_method_dir)
      args = ["bash", "process_1000g.sh"]
      print "Running process_1000g.sh with the following args: %s" % (args)
      print "It may take several minutes to download files... you will begin to see output once the downloads are complete."
      sp = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = sp.communicate()
      if out:
          print "standard output of subprocess:"
          print out
      if err:
          print "standard error of subprocess:"
          print err
      print "Completed 1000 Genomes download and extraction."

class DownloadAndExtractFilesFromEXAC(luigi.Task):

    def run(self):
      os.chdir(exac_file_dir)

      # Download vcf.gz file
      exac_vcf_gz_url = "ftp://ftp.broadinstitute.org/pub/ExAC_release/current/subsets/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz"
      exac_vcf_gz_file_name = exac_vcf_gz_url.split('/')[-1]
      u = urllib2.urlopen(exac_vcf_gz_url)
      f = open(exac_vcf_gz_file_name, 'wb')
      meta = u.info()
      file_size = int(meta.getheaders("Content-Length")[0])
      print "Downloading: %s Bytes: %s" % (exac_vcf_gz_file_name, file_size)

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
      print "Finished downloading %s" % (exac_vcf_gz_file_name)

      # Download vcf.gz.tbi file
      exac_vcf_gz_tbi_url = "ftp://ftp.broadinstitute.org/pub/ExAC_release/current/subsets/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz.tbi"
      exac_vcf_gz_tbi_file_name = exac_vcf_gz_tbi_url.split('/')[-1]
      u = urllib2.urlopen(exac_vcf_gz_tbi_url)
      f = open(exac_vcf_gz_tbi_file_name, 'wb')
      meta = u.info()
      file_size = int(meta.getheaders("Content-Length")[0])
      print "Downloading: %s Bytes: %s" % (exac_vcf_gz_tbi_file_name, file_size)

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
      print "Finished downloading %s" % (exac_vcf_gz_tbi_file_name)

      # tabix -h $EXAC/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz 17:41191488-41322420 > $EXAC/exac.brca1.hg19.vcf
      args = ["tabix", "-h", exac_file_dir + "/" + exac_vcf_gz_file_name, "17:41191488-41322420"]
      exac_brca1_hg19_vcf_file = exac_file_dir + "/exac.brca1.hg19.vcf"
      writable_exac_brca1_hg19_vcf_file = open(exac_brca1_hg19_vcf_file, 'w')
      print "Running tabix with the following args: %s" % (args)
      sp = subprocess.Popen(args, stdout=writable_exac_brca1_hg19_vcf_file, stderr=subprocess.PIPE)
      out, err = sp.communicate()
      if out:
          print "standard output of subprocess:"
          print out
      if err:
          print "standard error of subprocess:"
          print err
      print "Completed writing %s." % (exac_brca1_hg19_vcf_file)

      # tabix -h $EXAC/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz 13:32889080-32973809 > $EXAC/exac.brca2.hg19.vcf
      args = ["tabix", "-h", exac_file_dir + "/" + exac_vcf_gz_file_name, "13:32889080-32973809"]
      exac_brca2_hg19_vcf_file = exac_file_dir + "/exac.brca2.hg19.vcf"
      writable_exac_brca2_hg19_vcf_file = open(exac_brca2_hg19_vcf_file, 'w')
      print "Running tabix with the following args: %s" % (args)
      sp = subprocess.Popen(args, stdout=writable_exac_brca2_hg19_vcf_file, stderr=subprocess.PIPE)
      out, err = sp.communicate()
      if out:
          print "standard output of subprocess:"
          print out
      if err:
          print "standard error of subprocess:"
          print err
      print "Completed writing %s." % (exac_brca2_hg19_vcf_file)

      # vcf-concat $EXAC/exac.brca1.hg19.vcf $EXAC/exac.brca2.hg19.vcf > $EXAC/exac.brca12.hg19.vcf
      args = ["vcf-concat", exac_file_dir + "/exac.brca1.hg19.vcf", exac_file_dir + "/exac.brca2.hg19.vcf"]
      exac_brca12_hg19_vcf_file = exac_file_dir + "/exac.brca12.hg19.vcf"
      writable_exac_brca12_hg19_vcf_file = open(exac_brca12_hg19_vcf_file, 'w')
      print "Running tabix with the following args: %s" % (args)
      sp = subprocess.Popen(args, stdout=writable_exac_brca12_hg19_vcf_file, stderr=subprocess.PIPE)
      out, err = sp.communicate()
      if out:
          print "standard output of subprocess:"
          print out
      if err:
          print "standard error of subprocess:"
          print err
      print "Completed concatenation of exac brca1/2 data into %s." % (exac_brca12_hg19_vcf_file)

      # CrossMap.py vcf $BRCA_RESOURCES/hg19ToHg38.over.chain.gz $EXAC/exac.brca12.hg19.vcf $BRCA_RESOURCES/hg38.fa $EXAC/exac.brca12.hg38.vcf
      args = ["CrossMap.py", "vcf", brca_resources_dir + "/hg19ToHg38.over.chain.gz", exac_file_dir + "/exac.brca12.hg19.vcf", brca_resources_dir + "/hg38.fa", exac_file_dir + "/exac.brca12.hg38.vcf"]
      print "Running CrossMap.py with the following args: %s" % (args)
      sp = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = sp.communicate()
      if out:
          print "standard output of subprocess:"
          print out
      if err:
          print "standard error of subprocess:"
          print err
      print "Completed crossmapping hg19 and hg38."

      # vcf-sort $EXAC/exac.brca12.hg38.vcf > $EXAC/exac.brca12.sorted.hg38.vcf
      args = ["vcf-sort", exac_file_dir + "/exac.brca12.hg38.vcf"]
      exac_brca12_sorted_hg19_vcf_file = exac_file_dir + "/exac.brca12.sorted.hg38.vcf"
      writable_exac_brca12_sorted_hg19_vcf_file = open(exac_brca12_sorted_hg19_vcf_file, 'w')
      print "Running tabix with the following args: %s" % (args)
      sp = subprocess.Popen(args, stdout=writable_exac_brca12_sorted_hg19_vcf_file, stderr=subprocess.PIPE)
      out, err = sp.communicate()
      if out:
          print "standard output of subprocess:"
          print out
      if err:
          print "standard error of subprocess:"
          print err
      print "Completed sorting of exac data into %s." % (exac_brca12_sorted_hg19_vcf_file)

class RunAll(luigi.WrapperTask):
    date = luigi.DateParameter(default=datetime.date.today())
    u = luigi.Parameter()
    p = luigi.Parameter()
    def requires(self):
        # yield ConvertLatestClinvarToVCF(self.date)
        # yield DownloadAndExtractFilesFromESPTar(self.date)
        yield DownloadAndExtractFilesFromBIC(self.date, self.u, self.p)
        # yield DownloadAndExtractFilesFromG1K(self.date)
        # yield DownloadAndExtractFilesFromEXAC(self.date)
        # yield ExtractAndConvertFilesFromEXLOVD(self.date)
        # yield ExtractAndConvertFilesFromLOVD(self.date)

    def output(self):
        return luigi.LocalTarget('RunAll.txt')