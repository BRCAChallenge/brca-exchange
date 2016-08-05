import subprocess
import os
import urllib2
import luigi

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
