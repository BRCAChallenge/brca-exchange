import datetime
import json
import os
import shutil
import subprocess
import tarfile
import tempfile
from shutil import copy

import luigi
from luigi.util import requires

import esp_processing
import gnomad_processing
import pipeline_common
import pipeline_utils
from pipeline_common import DefaultPipelineTask

#######################################
# Default Globals / Env / Directories #
#######################################

luigi_dir = os.getcwd()

clinvar_method_dir = os.path.abspath('../clinvar')
lovd_method_dir = os.path.abspath('../lovd')
g1k_method_dir = os.path.abspath('../1000_Genomes')
enigma_method_dir = os.path.abspath('../enigma')
functional_assays_method_dir = os.path.abspath('../functional_assays')
data_merging_method_dir = os.path.abspath('../data_merging')
priors_method_dir = os.path.abspath('../splicing')
priors_filter_method_dir = os.path.abspath('../splicingfilter')
utilities_method_dir = os.path.abspath('../utilities')
vr_method_dir = os.path.abspath('../vr')


###############################################
#                   CLINVAR                   #
###############################################


class DownloadLatestClinvarData(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(
            self.cfg.file_parent_dir + "/ClinVar/ClinVarFullRelease_00-latest.xml.gz")

    def run(self):
        clinvar_file_dir = pipeline_utils.create_path_if_nonexistent(
            self.cfg.file_parent_dir + "/ClinVar")
        os.chdir(clinvar_file_dir)

        clinvar_data_url = "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz"
        pipeline_utils.download_file_and_display_progress(clinvar_data_url)


@requires(DownloadLatestClinvarData)
class ConvertLatestClinvarDataToXML(DefaultPipelineTask):

    def output(self):
        artifacts_dir = pipeline_utils.create_path_if_nonexistent(
            self.cfg.output_dir + "/release/artifacts/")
        return luigi.LocalTarget(
            self.cfg.file_parent_dir + "/ClinVar/ClinVarBrca.xml")

    def run(self):
        artifacts_dir = pipeline_utils.create_path_if_nonexistent(
            self.cfg.output_dir + "/release/artifacts/")
        clinvar_file_dir = self.cfg.file_parent_dir + "/ClinVar"
        os.chdir(clinvar_method_dir)

        clinvar_xml_file = clinvar_file_dir + "/ClinVarBrca.xml"
        writable_clinvar_xml_file = open(clinvar_xml_file, "w")
        args = ["python", "filter_clinvar_brca.py", self.input().path,
                self.output().path]
        print "Running clinVarBrca.py with the following args: %s. This takes a while..." % (
            args)
        sp = subprocess.Popen(args, stdout=writable_clinvar_xml_file,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(clinvar_xml_file)


@requires(ConvertLatestClinvarDataToXML)
class ConvertClinvarXMLToTXT(DefaultPipelineTask):

    def output(self):
        return luigi.LocalTarget(
            self.cfg.file_parent_dir + "/ClinVar/ClinVarBrca.txt")

    def run(self):
        clinvar_file_dir = self.cfg.file_parent_dir + "/ClinVar"
        os.chdir(clinvar_method_dir)

        clinvar_txt_file = clinvar_file_dir + "/ClinVarBrca.txt"
        writable_clinvar_txt_file = open(clinvar_txt_file, "w")
        args = ["python", "clinVarParse.py",
                clinvar_file_dir + "/ClinVarBrca.xml",
                "--logs", clinvar_file_dir + "/clinvar_xml_to_txt.log",
                "--assembly", "GRCh38"]
        print "Running clinVarParse.py with the following args: %s" % (args)
        sp = subprocess.Popen(args, stdout=writable_clinvar_txt_file,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(clinvar_txt_file)


@requires(ConvertClinvarXMLToTXT)
class ConvertClinvarTXTToVCF(DefaultPipelineTask):

    def output(self):
        return luigi.LocalTarget(
            self.cfg.file_parent_dir + "/ClinVar/ClinVarBrca.vcf")

    def run(self):
        clinvar_file_dir = self.cfg.file_parent_dir + "/ClinVar"
        clinvar_vcf_file = clinvar_file_dir + "/ClinVarBrca.vcf"

        os.chdir(data_merging_method_dir)
        args = ["python", "convert_tsv_to_vcf.py", "-i",
                clinvar_file_dir + "/ClinVarBrca.txt", "-o",
                clinvar_file_dir + "/ClinVarBrca.vcf", "-s", "ClinVar"]
        print "Running convert_tsv_to_vcf.py with the following args: %s" % (
            args)
        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(clinvar_vcf_file)


@requires(ConvertClinvarTXTToVCF)
class CopyClinvarVCFToOutputDir(DefaultPipelineTask):

    def output(self):
        return luigi.LocalTarget(self.cfg.output_dir + "/ClinVarBrca.vcf")

    def run(self):
        clinvar_file_dir = self.cfg.file_parent_dir + "/ClinVar"
        pipeline_utils.create_path_if_nonexistent(self.cfg.output_dir)

        copy(self.cfg.file_parent_dir + "/ClinVar/ClinVarBrca.vcf",
             self.cfg.output_dir)
        pipeline_utils.check_file_for_contents(
            self.cfg.output_dir + "/ClinVarBrca.vcf")


###############################################
#                     BIC                     #
###############################################


class DownloadBICData(DefaultPipelineTask):
    def output(self):
        bic_file_dir = self.cfg.file_parent_dir + '/BIC'
        return luigi.LocalTarget(bic_file_dir + "/bic_brca12.sorted.hg38.vcf")

    def run(self):
        bic_file_dir = pipeline_utils.create_path_if_nonexistent(
            self.cfg.file_parent_dir + '/BIC')

        os.chdir(bic_file_dir)

        brca1_data_url = "https://brcaexchange.org/backend/downloads/bic_brca12.sorted.hg38.vcf"
        pipeline_utils.download_file_and_display_progress(brca1_data_url)


@requires(DownloadBICData)
class CopyBICOutputToOutputDir(DefaultPipelineTask):

    def output(self):
        bic_file_dir = self.cfg.file_parent_dir + '/BIC'
        return luigi.LocalTarget(
            self.cfg.output_dir + "/bic_brca12.sorted.hg38.vcf")

    def run(self):
        bic_file_dir = self.cfg.file_parent_dir + '/BIC'
        pipeline_utils.create_path_if_nonexistent(self.cfg.output_dir)

        copy(bic_file_dir + "/bic_brca12.sorted.hg38.vcf", self.cfg.output_dir)
        pipeline_utils.check_file_for_contents(
            self.cfg.output_dir + "/bic_brca12.sorted.hg38.vcf")


###############################################
#                  exLOVD                     #
###############################################


class ExtractDataFromLatestEXLOVD(DefaultPipelineTask):
    dir_name = 'exLOVD'

    def __init__(self):
        super(ExtractDataFromLatestEXLOVD, self).__init__()
        self.ex_lovd_file_dir = os.path.join(self.cfg.file_parent_dir,
                                        ExtractDataFromLatestEXLOVD.dir_name)

    def output(self):
        return {'brca1': luigi.LocalTarget(os.path.join(self.ex_lovd_file_dir, "BRCA1.txt")),
                'brca2': luigi.LocalTarget(os.path.join(self.ex_lovd_file_dir, "BRCA2.txt"))}

    def run(self):
        pipeline_utils.create_path_if_nonexistent(self.ex_lovd_file_dir)

        # calculating host path because we are running a docker within a docker
        ex_lovd_file_dir_host = os.path.join(os.path.dirname(self.cfg.output_dir_host), ExtractDataFromLatestEXLOVD.dir_name)

        os.chdir(lovd_method_dir)

        ex_lovd_data_host_url = "http://hci-exlovd.hci.utah.edu/"

        args = ['bash', 'extract_latest_exlovd.sh', ex_lovd_file_dir_host, "-u", ex_lovd_data_host_url, "-l", "BRCA1",
                "BRCA2", "-o", "/data"]

        print "Running extract_data.py with the following args: %s" % (args)
        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)
        print "Extracted data from %s." % (ex_lovd_data_host_url)


@requires(ExtractDataFromLatestEXLOVD)
class ConvertEXLOVDBRCA1ExtractToVCF(DefaultPipelineTask):

    def output(self):
        ex_lovd_file_dir = self.cfg.file_parent_dir + "/exLOVD"
        return luigi.LocalTarget(ex_lovd_file_dir + "/exLOVD_brca1.hg19.vcf")

    def run(self):
        ex_lovd_file_dir = self.cfg.file_parent_dir + "/exLOVD"
        brca_resources_dir = self.cfg.resources_dir
        artifacts_dir = pipeline_utils.create_path_if_nonexistent(
            self.cfg.output_dir + "/release/artifacts/")

        os.chdir(lovd_method_dir)

        args = ["./lovd2vcf.py", "-i", ex_lovd_file_dir + "/BRCA1.txt", "-o",
                ex_lovd_file_dir + "/exLOVD_brca1.hg19.vcf", "-a",
                "exLOVDAnnotation",
                "-r", brca_resources_dir + "/refseq_annotation.hg19.gp", "-g",
                brca_resources_dir + "/hg19.fa", "-e",
                artifacts_dir + "exLOVD_BRCA1_error_variants.txt",
                "-s", "exLOVD"]
        print "Running lovd2vcf with the following args: %s" % (args)
        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(
            ex_lovd_file_dir + "/exLOVD_brca1.hg19.vcf")


@requires(ConvertEXLOVDBRCA1ExtractToVCF)
class ConvertEXLOVDBRCA2ExtractToVCF(DefaultPipelineTask):

    def output(self):
        ex_lovd_file_dir = self.cfg.file_parent_dir + "/exLOVD"
        return luigi.LocalTarget(ex_lovd_file_dir + "/exLOVD_brca2.hg19.vcf")

    def run(self):
        ex_lovd_file_dir = self.cfg.file_parent_dir + "/exLOVD"
        brca_resources_dir = self.cfg.resources_dir
        artifacts_dir = pipeline_utils.create_path_if_nonexistent(
            self.cfg.output_dir + "/release/artifacts/")

        os.chdir(lovd_method_dir)

        args = ["./lovd2vcf.py", "-i", ex_lovd_file_dir + "/BRCA2.txt", "-o",
                ex_lovd_file_dir + "/exLOVD_brca2.hg19.vcf", "-a",
                "exLOVDAnnotation",
                "-r", brca_resources_dir + "/refseq_annotation.hg19.gp", "-g",
                brca_resources_dir + "/hg19.fa", "-e",
                artifacts_dir + "exLOVD_BRCA2_error_variants.txt",
                "-s", "exLOVD"]
        print "Running lovd2vcf with the following args: %s" % (args)
        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(
            ex_lovd_file_dir + "/exLOVD_brca2.hg19.vcf")


@requires(ConvertEXLOVDBRCA2ExtractToVCF)
class ConcatenateEXLOVDVCFFiles(DefaultPipelineTask):

    def output(self):
        ex_lovd_file_dir = self.cfg.file_parent_dir + "/exLOVD"
        return luigi.LocalTarget(ex_lovd_file_dir + "/exLOVD_brca12.hg19.vcf")

    def run(self):
        ex_lovd_file_dir = self.cfg.file_parent_dir + "/exLOVD"

        ex_lovd_brca12_hg19_vcf_file = ex_lovd_file_dir + "/exLOVD_brca12.hg19.vcf"
        writable_ex_lovd_brca12_hg19_vcf_file = open(
            ex_lovd_brca12_hg19_vcf_file, 'w')
        args = ["vcf-concat", ex_lovd_file_dir + "/exLOVD_brca1.hg19.vcf",
                ex_lovd_file_dir + "/exLOVD_brca2.hg19.vcf"]
        print "Running vcf-concat with the following args: %s" % (args)
        sp = subprocess.Popen(args,
                              stdout=writable_ex_lovd_brca12_hg19_vcf_file,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(
            ex_lovd_file_dir + "/exLOVD_brca12.hg19.vcf")


@requires(ConcatenateEXLOVDVCFFiles)
class CrossmapConcatenatedEXLOVDData(DefaultPipelineTask):

    def output(self):
        ex_lovd_file_dir = self.cfg.file_parent_dir + "/exLOVD"
        return luigi.LocalTarget(ex_lovd_file_dir + "/exLOVD_brca12.hg38.vcf")

    def run(self):
        ex_lovd_file_dir = self.cfg.file_parent_dir + "/exLOVD"
        brca_resources_dir = self.cfg.resources_dir

        args = ["CrossMap.py", "vcf",
                brca_resources_dir + "/hg19ToHg38.over.chain.gz",
                ex_lovd_file_dir + "/exLOVD_brca12.hg19.vcf",
                brca_resources_dir + "/hg38.fa",
                ex_lovd_file_dir + "/exLOVD_brca12.hg38.vcf"]
        print "Running CrossMap.py with the following args: %s" % (args)
        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(
            ex_lovd_file_dir + "/exLOVD_brca12.hg38.vcf")


@requires(CrossmapConcatenatedEXLOVDData)
class SortEXLOVDOutput(DefaultPipelineTask):

    def output(self):
        ex_lovd_file_dir = self.cfg.file_parent_dir + "/exLOVD"
        return luigi.LocalTarget(
            ex_lovd_file_dir + "/exLOVD_brca12.sorted.hg38.vcf")

    def run(self):
        ex_lovd_file_dir = self.cfg.file_parent_dir + "/exLOVD"

        sorted_ex_lovd_output_file = ex_lovd_file_dir + "/exLOVD_brca12.sorted.hg38.vcf"
        writable_sorted_ex_lovd_output_file = open(sorted_ex_lovd_output_file,
                                                   'w')
        args = ["vcf-sort", ex_lovd_file_dir + "/exLOVD_brca12.hg38.vcf"]
        print "Running vcf-sort with the following args: %s" % (args)
        sp = subprocess.Popen(args, stdout=writable_sorted_ex_lovd_output_file,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)
        print "Sorted BRCA1/2 hg38 vcf file into %s" % (
            writable_sorted_ex_lovd_output_file)

        pipeline_utils.check_file_for_contents(
            ex_lovd_file_dir + "/exLOVD_brca12.sorted.hg38.vcf")


@requires(SortEXLOVDOutput)
class CopyEXLOVDOutputToOutputDir(DefaultPipelineTask):

    def output(self):
        return luigi.LocalTarget(
            self.cfg.output_dir + "/exLOVD_brca12.sorted.hg38.vcf")

    def run(self):
        ex_lovd_file_dir = self.cfg.file_parent_dir + "/exLOVD"
        pipeline_utils.create_path_if_nonexistent(self.cfg.output_dir)

        copy(ex_lovd_file_dir + "/exLOVD_brca12.sorted.hg38.vcf",
             self.cfg.output_dir)

        pipeline_utils.check_file_for_contents(
            self.cfg.output_dir + "/exLOVD_brca12.sorted.hg38.vcf")


##############################################
#               sharedLOVD                   #
##############################################


class DownloadLOVDInputFile(DefaultPipelineTask):
    """ Downloads the shared LOVD data

    If the pipeline is run on a machine from which it is not possible to download the data (currently IP based authentication)
    the file can be manually staged in the path of `lovd_data_file`. In this case, the task will not be run.
    """

    lovd_data_file = luigi.Parameter(default='',
                                     description='path, where the shared LOVD data will be stored')

    shared_lovd_data_url = luigi.Parameter(
        default='https://databases.lovd.nl/shared/export/BRCA',
        description='URL to download shared LOVD data from')

    def output(self):
        if len(str(self.lovd_data_file)) == 0:
            path = self.cfg.file_parent_dir + "/LOVD/BRCA.txt"
        else:
            path = str(self.lovd_data_file)

        return luigi.LocalTarget(path)

    def run(self):
        pipeline_utils.create_path_if_nonexistent(
            os.path.dirname(self.output().path))
        data = pipeline_utils.urlopen_with_retry(
            self.shared_lovd_data_url).read()
        with open(self.output().path, "wb") as f:
            f.write(data)


@requires(DownloadLOVDInputFile)
class NormalizeLOVDSubmissions(DefaultPipelineTask):

    def output(self):
        return luigi.LocalTarget(
            self.cfg.file_parent_dir + "/LOVD/LOVD_normalized.tsv")

    def run(self):
        artifacts_dir = pipeline_utils.create_path_if_nonexistent(
            self.cfg.output_dir + "/release/artifacts")

        os.chdir(lovd_method_dir)

        args = ["python", "normalizeLOVDSubmissions.py", "-i",
                self.input().path, "-o",
                self.output().path]

        print "Running NormalizeLOVDSubmissions with the following args: %s" % (args)

        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(self.output().path)


@requires(NormalizeLOVDSubmissions)
class CombineEquivalentLOVDSubmissions(DefaultPipelineTask):

    def output(self):
        return luigi.LocalTarget(
            self.cfg.file_parent_dir + "/LOVD/LOVD_normalized_combined.tsv")

    def run(self):
        artifacts_dir = pipeline_utils.create_path_if_nonexistent(
            self.cfg.output_dir + "/release/artifacts")

        os.chdir(lovd_method_dir)

        args = ["python", "combineEquivalentVariantSubmissions.py", "-i",
                self.input().path, "-o",
                self.output().path]

        print "Running combineEquivalentVariantSubmissions.py with the following args: %s" % (
            args)

        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(self.output().path)


@requires(CombineEquivalentLOVDSubmissions)
class ConvertSharedLOVDToVCF(DefaultPipelineTask):

    def output(self):
        return luigi.LocalTarget(
            self.cfg.file_parent_dir + "/LOVD/sharedLOVD_brca12.hg19.vcf")

    def run(self):
        brca_resources_dir = self.cfg.resources_dir
        artifacts_dir = pipeline_utils.create_path_if_nonexistent(
            self.cfg.output_dir + "/release/artifacts")

        os.chdir(lovd_method_dir)

        args = ["python", "lovd2vcf.py", "-i", self.input().path, "-o",
                self.output().path, "-a", "sharedLOVDAnnotation",
                "-r", brca_resources_dir + "/refseq_annotation.hg19.gp", "-g",
                brca_resources_dir + "/hg19.fa", "-e",
                artifacts_dir + "/LOVD_error_variants.txt",
                "-s", "LOVD"]

        print "Running lovd2vcf with the following args: %s" % (args)

        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(self.output().path)


@requires(ConvertSharedLOVDToVCF)
class CrossmapConcatenatedSharedLOVDData(DefaultPipelineTask):

    def output(self):
        lovd_file_dir = self.cfg.file_parent_dir + "/LOVD"
        return luigi.LocalTarget(lovd_file_dir + "/sharedLOVD_brca12.hg38.vcf")

    def run(self):
        lovd_file_dir = self.cfg.file_parent_dir + "/LOVD"
        brca_resources_dir = self.cfg.resources_dir

        args = ["CrossMap.py", "vcf",
                brca_resources_dir + "/hg19ToHg38.over.chain.gz",
                lovd_file_dir + "/sharedLOVD_brca12.hg19.vcf",
                brca_resources_dir + "/hg38.fa",
                lovd_file_dir + "/sharedLOVD_brca12.hg38.vcf"]
        print "Running CrossMap.py with the following args: %s" % (args)
        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(
            lovd_file_dir + "/sharedLOVD_brca12.hg38.vcf")


@requires(CrossmapConcatenatedSharedLOVDData)
class SortSharedLOVDOutput(DefaultPipelineTask):

    def output(self):
        lovd_file_dir = self.cfg.file_parent_dir + "/LOVD"
        return luigi.LocalTarget(
            lovd_file_dir + "/sharedLOVD_brca12.sorted.hg38.vcf")

    def run(self):
        lovd_file_dir = self.cfg.file_parent_dir + "/LOVD"

        sorted_lovd_output_file = lovd_file_dir + "/sharedLOVD_brca12.sorted.hg38.vcf"
        writable_sorted_lovd_output_file = open(sorted_lovd_output_file, 'w')
        args = ["vcf-sort", lovd_file_dir + "/sharedLOVD_brca12.hg38.vcf"]
        print "Running vcf-sort with the following args: %s" % (args)
        sp = subprocess.Popen(args, stdout=writable_sorted_lovd_output_file,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)
        print "Sorted BRCA1/2 hg38 vcf file into %s" % (
            writable_sorted_lovd_output_file)

        pipeline_utils.check_file_for_contents(
            lovd_file_dir + "/sharedLOVD_brca12.sorted.hg38.vcf")


@requires(SortSharedLOVDOutput)
class CopySharedLOVDOutputToOutputDir(DefaultPipelineTask):

    def output(self):
        return luigi.LocalTarget(
            self.cfg.output_dir + "/sharedLOVD_brca12.sorted.hg38.vcf")

    def run(self):
        lovd_file_dir = self.cfg.file_parent_dir + "/LOVD"
        pipeline_utils.create_path_if_nonexistent(self.cfg.output_dir)

        copy(lovd_file_dir + "/sharedLOVD_brca12.sorted.hg38.vcf",
             self.cfg.output_dir)

        pipeline_utils.check_file_for_contents(
            self.cfg.output_dir + "/sharedLOVD_brca12.sorted.hg38.vcf")


###############################################
#                    G1K                      #
###############################################


class DownloadG1KCHR13GZ(DefaultPipelineTask):
    def output(self):
        g1k_file_dir = self.cfg.file_parent_dir + '/G1K'
        return luigi.LocalTarget(
            g1k_file_dir + "/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")

    def run(self):
        g1k_file_dir = pipeline_utils.create_path_if_nonexistent(
            self.cfg.file_parent_dir + '/G1K')

        os.chdir(g1k_file_dir)

        chr13_vcf_gz_url = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
        pipeline_utils.download_file_and_display_progress(chr13_vcf_gz_url)


@requires(DownloadG1KCHR13GZ)
class DownloadG1KCHR17GZ(DefaultPipelineTask):

    def output(self):
        g1k_file_dir = self.cfg.file_parent_dir + '/G1K'
        return luigi.LocalTarget(
            g1k_file_dir + "/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")

    def run(self):
        g1k_file_dir = self.cfg.file_parent_dir + '/G1K'
        os.chdir(g1k_file_dir)

        chr17_vcf_gz_url = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
        pipeline_utils.download_file_and_display_progress(chr17_vcf_gz_url)


@requires(DownloadG1KCHR17GZ)
class DownloadG1KCHR13GZTBI(DefaultPipelineTask):

    def output(self):
        g1k_file_dir = self.cfg.file_parent_dir + '/G1K'
        return luigi.LocalTarget(
            g1k_file_dir + "/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi")

    def run(self):
        g1k_file_dir = self.cfg.file_parent_dir + '/G1K'
        os.chdir(g1k_file_dir)

        chr13_vcf_gz_tbi_url = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi"
        pipeline_utils.download_file_and_display_progress(chr13_vcf_gz_tbi_url)


@requires(DownloadG1KCHR13GZTBI)
class DownloadG1KCHR17GZTBI(DefaultPipelineTask):

    def output(self):
        g1k_file_dir = self.cfg.file_parent_dir + '/G1K'
        return luigi.LocalTarget(
            g1k_file_dir + "/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi")

    def run(self):
        g1k_file_dir = self.cfg.file_parent_dir + '/G1K'
        os.chdir(g1k_file_dir)

        chr17_vcf_gz_tbi_url = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi"
        pipeline_utils.download_file_and_display_progress(chr17_vcf_gz_tbi_url)


@requires(DownloadG1KCHR17GZTBI)
class ExtractCHR13BRCAData(DefaultPipelineTask):

    def output(self):
        g1k_file_dir = self.cfg.file_parent_dir + '/G1K'
        return luigi.LocalTarget(g1k_file_dir + "/chr13_brca2_1000g_GRCh37.vcf")

    def run(self):
        g1k_file_dir = self.cfg.file_parent_dir + '/G1K'

        chr13_brca2_vcf_file = g1k_file_dir + "/chr13_brca2_1000g_GRCh37.vcf"
        writable_chr13_brca2_vcf_file = open(chr13_brca2_vcf_file, "w")
        args = ["tabix", "-h",
                g1k_file_dir + "/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
                "13:32889617-32973809"]
        print "Running tabix with the following args: %s" % (args)
        sp = subprocess.Popen(args, stdout=writable_chr13_brca2_vcf_file,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(chr13_brca2_vcf_file)


@requires(ExtractCHR13BRCAData)
class ExtractCHR17BRCAData(DefaultPipelineTask):

    def output(self):
        g1k_file_dir = self.cfg.file_parent_dir + '/G1K'
        return luigi.LocalTarget(g1k_file_dir + "/chr17_brca1_1000g_GRCh37.vcf")

    def run(self):
        g1k_file_dir = self.cfg.file_parent_dir + '/G1K'

        chr17_brca1_vcf_file = g1k_file_dir + "/chr17_brca1_1000g_GRCh37.vcf"
        writable_chr17_brca1_vcf_file = open(chr17_brca1_vcf_file, "w")
        args = ["tabix", "-h",
                g1k_file_dir + "/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
                "17:41196312-41277500"]
        print "Running tabix with the following args: %s" % (args)
        sp = subprocess.Popen(args, stdout=writable_chr17_brca1_vcf_file,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(chr17_brca1_vcf_file)


@requires(ExtractCHR17BRCAData)
class ConcatenateG1KData(DefaultPipelineTask):

    def output(self):
        g1k_file_dir = self.cfg.file_parent_dir + '/G1K'
        return luigi.LocalTarget(g1k_file_dir + "/brca12_1000g_GRCh37.vcf")

    def run(self):
        g1k_file_dir = self.cfg.file_parent_dir + '/G1K'
        concatenated_g1k_vcf = g1k_file_dir + "/brca12_1000g_GRCh37.vcf"
        writable_concatenated_g1k_vcf = open(concatenated_g1k_vcf, "w")
        args = ["vcf-concat", g1k_file_dir + "/chr13_brca2_1000g_GRCh37.vcf",
                g1k_file_dir + "/chr17_brca1_1000g_GRCh37.vcf"]
        print "Running vcf-concat with the following args: %s" % (args)
        sp = subprocess.Popen(args, stdout=writable_concatenated_g1k_vcf,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(concatenated_g1k_vcf)


@requires(ConcatenateG1KData)
class CrossmapConcatenatedG1KData(DefaultPipelineTask):

    def output(self):
        g1k_file_dir = self.cfg.file_parent_dir + '/G1K'
        return luigi.LocalTarget(g1k_file_dir + "/1000G_brca.hg38.vcf")

    def run(self):
        g1k_file_dir = self.cfg.file_parent_dir + '/G1K'
        brca_resources_dir = self.cfg.resources_dir

        args = ["CrossMap.py", "vcf",
                brca_resources_dir + "/hg19ToHg38.over.chain.gz",
                g1k_file_dir + "/brca12_1000g_GRCh37.vcf",
                brca_resources_dir + "/hg38.fa",
                g1k_file_dir + "/1000G_brca.hg38.vcf"]
        print "Running CrossMap.py with the following args: %s" % (args)
        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(
            g1k_file_dir + "/1000G_brca.hg38.vcf")


@requires(CrossmapConcatenatedG1KData)
class SortG1KData(DefaultPipelineTask):

    def output(self):
        g1k_file_dir = self.cfg.file_parent_dir + '/G1K'
        return luigi.LocalTarget(g1k_file_dir + "/1000G_brca.sorted.hg38.vcf")

    def run(self):
        g1k_file_dir = self.cfg.file_parent_dir + '/G1K'

        sorted_g1k_output_file = g1k_file_dir + "/1000G_brca.sorted.hg38.vcf"
        writable_sorted_g1k_output_file = open(sorted_g1k_output_file, 'w')
        args = ["vcf-sort", g1k_file_dir + "/1000G_brca.hg38.vcf"]
        print "Running vcf-sort with the following args: %s" % (args)
        sp = subprocess.Popen(args, stdout=writable_sorted_g1k_output_file,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(
            g1k_file_dir + "/1000G_brca.sorted.hg38.vcf")


@requires(SortG1KData)
class CopyG1KOutputToOutputDir(DefaultPipelineTask):

    def output(self):
        return luigi.LocalTarget(
            self.cfg.output_dir + "/1000G_brca.sorted.hg38.vcf")

    def run(self):
        g1k_file_dir = self.cfg.file_parent_dir + '/G1K'
        pipeline_utils.create_path_if_nonexistent(self.cfg.output_dir)

        copy(g1k_file_dir + "/1000G_brca.sorted.hg38.vcf", self.cfg.output_dir)

        pipeline_utils.check_file_for_contents(
            self.cfg.output_dir + "/1000G_brca.sorted.hg38.vcf")


###############################################
#                    EXAC                     #
###############################################


class DownloadStaticExACData(DefaultPipelineTask):
    def output(self):
        exac_file_dir = self.cfg.file_parent_dir + '/exac'
        return luigi.LocalTarget(
            exac_file_dir + "/exac.brca12.sorted.hg38.vcf")

    def run(self):
        exac_file_dir = pipeline_utils.create_path_if_nonexistent(
            self.cfg.file_parent_dir + '/exac')

        os.chdir(exac_file_dir)

        exac_vcf_gz_url = "https://brcaexchange.org/backend/downloads/exac.brca12.sorted.hg38.vcf"
        pipeline_utils.download_file_and_display_progress(exac_vcf_gz_url)


@requires(DownloadStaticExACData)
class CopyEXACOutputToOutputDir(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(
            self.cfg.output_dir + "/exac.brca12.sorted.hg38.vcf")

    def run(self):
        exac_file_dir = self.cfg.file_parent_dir + '/exac'
        pipeline_utils.create_path_if_nonexistent(self.cfg.output_dir)

        copy(exac_file_dir + "/exac.brca12.sorted.hg38.vcf",
             self.cfg.output_dir)

        pipeline_utils.check_file_for_contents(
            self.cfg.output_dir + "/exac.brca12.sorted.hg38.vcf")


###############################################
#                  ENIGMA                     #
###############################################

@requires(ConvertLatestClinvarDataToXML)
class FilterEnigmaAssertions(DefaultPipelineTask):
    def output(self):
        out_path = os.path.join(self.cfg.file_parent_dir, 'enigma',
                                'enigma_clinvar.xml')
        return luigi.LocalTarget(out_path)

    def run(self):
        pipeline_utils.create_path_if_nonexistent(
            os.path.join(self.cfg.file_parent_dir, 'enigma'))
        os.chdir(clinvar_method_dir)

        args = ["python", "filter_enigma_data.py", self.input().path,
                self.output().path]

        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)


@requires(FilterEnigmaAssertions)
class ExtractEnigmaFromClinvar(DefaultPipelineTask):
    def output(self):
        out_path = os.path.join(self.cfg.file_parent_dir, 'enigma',
                                'enigma_from_clinvar.tsv')
        return luigi.LocalTarget(out_path)

    def run(self):
        os.chdir(clinvar_method_dir)

        args = ["python", "enigma_from_clinvar.py", self.input().path,
                self.output().path,
                '--logs', os.path.join(self.cfg.file_parent_dir, 'enigma', 'enigma_from_clinvar.log')
               ]

        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        copy(self.output().path, self.cfg.output_dir)


@requires(ExtractEnigmaFromClinvar)
class CopyEnigmaOutputToOutputDir(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(self.cfg.output_dir + "/enigma_from_clinvar.tsv")

    def run(self):
        enigma_file_dir = self.cfg.file_parent_dir + '/enigma'
        pipeline_utils.create_path_if_nonexistent(self.cfg.output_dir)

        copy(enigma_file_dir + "/enigma_from_clinvar.tsv", self.cfg.output_dir)

        pipeline_utils.check_file_for_contents(self.cfg.output_dir + "/enigma_from_clinvar.tsv")



###############################################
#             FUNCTIONAL ASSAYS               #
###############################################


class DownloadFindlayBRCA1RingFunctionScoresInputFile(DefaultPipelineTask):
    findlay_BRCA1_ring_function_scores_url = luigi.Parameter(default='https://brcaexchange.org/backend/downloads/findlay_BRCA1_ring_function_scores.tsv',
                                            description='URL to download findlay_BRCA1_ring_function_scores data from')

    def output(self):
        return luigi.LocalTarget(self.cfg.file_parent_dir + "/functional_assays/findlay_BRCA1_ring_function_scores.tsv")

    def run(self):
        pipeline_utils.create_path_if_nonexistent(os.path.dirname(self.output().path))
        data = pipeline_utils.urlopen_with_retry(self.findlay_BRCA1_ring_function_scores_url).read()
        with open(self.output().path, "wb") as f:
            f.write(data)


@requires(DownloadFindlayBRCA1RingFunctionScoresInputFile)
class ParseFindlayBRCA1RingFunctionScores(DefaultPipelineTask):

    def output(self):
        return luigi.LocalTarget(self.cfg.file_parent_dir + "/functional_assays/findlay_BRCA1_ring_function_scores.clean.tsv")

    def run(self):
        os.chdir(functional_assays_method_dir)

        args = ["python", "parse_functional_assay_data.py", "-i", self.input().path, "-o",
                self.output().path]

        print "Running parse_functional_assay_data.py with the following args: %s" % (args)

        sp = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(self.output().path)


@requires(ParseFindlayBRCA1RingFunctionScores)
class ConvertFindlayBRCA1RingFunctionScoresToVCF(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(self.cfg.file_parent_dir + "/functional_assays/findlay_BRCA1_ring_function_scores.clean.hg19.vcf")

    def run(self):
        brca_resources_dir = self.cfg.resources_dir
        artifacts_dir = pipeline_utils.create_path_if_nonexistent(self.cfg.output_dir + "/release/artifacts")

        os.chdir(functional_assays_method_dir)

        args = ["python", "functional_assays_to_vcf.py", "-i", self.input().path, "-o",
                self.output().path, "-a", "functionalAssayAnnotation",
                "-r", brca_resources_dir + "/refseq_annotation.hg19.gp", "-g",
                brca_resources_dir + "/hg19.fa", "-l", artifacts_dir + "/findlay_BRCA1_ring_function_scores_error_variants.log",
                "-s", "FindlayBRCA1RingFunctionScores"]

        print "Running functional_assays_to_vcf with the following args: %s" % (args)

        sp = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(self.output().path)


@requires(ConvertFindlayBRCA1RingFunctionScoresToVCF)
class CrossmapFindlayBRCA1RingFunctionScores(DefaultPipelineTask):

    def output(self):
        functional_assays_file_dir = self.cfg.file_parent_dir + "/functional_assays"
        return luigi.LocalTarget(functional_assays_file_dir + "/findlay_BRCA1_ring_function_scores.clean.hg38.vcf")

    def run(self):
        functional_assays_file_dir = self.cfg.file_parent_dir + "/functional_assays"
        brca_resources_dir = self.cfg.resources_dir

        args = ["CrossMap.py", "vcf", brca_resources_dir + "/hg19ToHg38.over.chain.gz",
                functional_assays_file_dir + "/findlay_BRCA1_ring_function_scores.clean.hg19.vcf", brca_resources_dir + "/hg38.fa",
                functional_assays_file_dir + "/findlay_BRCA1_ring_function_scores.clean.hg38.vcf"]
        print "Running CrossMap.py with the following args: %s" % (args)
        sp = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(functional_assays_file_dir + "/findlay_BRCA1_ring_function_scores.clean.hg38.vcf")


@requires(CrossmapFindlayBRCA1RingFunctionScores)
class SortFindlayBRCA1RingFunctionScores(DefaultPipelineTask):

    def output(self):
        functional_assays_file_dir = self.cfg.file_parent_dir + "/functional_assays"
        return luigi.LocalTarget(functional_assays_file_dir + "/findlay_BRCA1_ring_function_scores.clean.sorted.hg38.vcf")

    def run(self):
        functional_assays_file_dir = self.cfg.file_parent_dir + "/functional_assays"

        sorted_findlay_BRCA1_ring_function_scores = functional_assays_file_dir + "/findlay_BRCA1_ring_function_scores.clean.sorted.hg38.vcf"
        writable_sorted_findlay_BRCA1_ring_function_scores = open(sorted_findlay_BRCA1_ring_function_scores, 'w')
        args = ["vcf-sort", functional_assays_file_dir + "/findlay_BRCA1_ring_function_scores.clean.hg38.vcf"]
        print "Running vcf-sort with the following args: %s" % (args)
        sp = subprocess.Popen(args, stdout=writable_sorted_findlay_BRCA1_ring_function_scores, stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)
        print "Sorted hg38 vcf file into %s" % (writable_sorted_findlay_BRCA1_ring_function_scores)

        pipeline_utils.check_file_for_contents(functional_assays_file_dir + "/findlay_BRCA1_ring_function_scores.clean.sorted.hg38.vcf")


@requires(SortFindlayBRCA1RingFunctionScores)
class CopyFindlayBRCA1RingFunctionScoresOutputToOutputDir(DefaultPipelineTask):

    def output(self):
        return luigi.LocalTarget(self.cfg.output_dir + "/findlay_BRCA1_ring_function_scores.clean.sorted.hg38.vcf")

    def run(self):
        functional_assays_file_dir = self.cfg.file_parent_dir + "/functional_assays"
        pipeline_utils.create_path_if_nonexistent(self.cfg.output_dir)

        copy(functional_assays_file_dir + "/findlay_BRCA1_ring_function_scores.clean.sorted.hg38.vcf", self.cfg.output_dir)

        pipeline_utils.check_file_for_contents(self.cfg.output_dir + "/findlay_BRCA1_ring_function_scores.clean.sorted.hg38.vcf")


###############################################
#            VARIANT COMPILATION              #
###############################################


class MergeVCFsIntoTSVFile(DefaultPipelineTask):
    def requires(self):
        yield pipeline_common.CopyOutputToOutputDir(self.cfg.output_dir,
                                                    esp_processing.SortConcatenatedESPData())
        yield pipeline_common.CopyOutputToOutputDir(self.cfg.output_dir,
                                                    gnomad_processing.SortGnomADData())
        yield CopyClinvarVCFToOutputDir()
        yield CopyBICOutputToOutputDir()
        yield CopyG1KOutputToOutputDir()
        yield CopyEXACOutputToOutputDir()
        yield CopyEXLOVDOutputToOutputDir()
        yield CopySharedLOVDOutputToOutputDir()
        yield CopyEnigmaOutputToOutputDir()
        yield CopyFindlayBRCA1RingFunctionScoresOutputToOutputDir()

    def output(self):
        artifacts_dir = pipeline_utils.create_path_if_nonexistent(
            self.cfg.output_dir + "/release/artifacts/")
        return [luigi.LocalTarget(artifacts_dir + "merged.tsv"),
                luigi.LocalTarget(artifacts_dir + "reports.tsv")]

    def run(self):
        print("running merged")
        print(self.input())
        artifacts_dir = pipeline_utils.create_path_if_nonexistent(
            self.cfg.output_dir + "/release/artifacts/")

        os.chdir(data_merging_method_dir)

        args = ["python", "variant_merging.py", "-i", self.cfg.output_dir + "/",
                "-o",
                artifacts_dir, "-a",
                artifacts_dir, "-v",
                "-c", self.cfg.gene_config_path]
        print "Running variant_merging.py with the following args: %s" % (args)
        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(artifacts_dir + "merged.tsv")
        pipeline_utils.check_file_for_contents(artifacts_dir + "reports.tsv")


@requires(MergeVCFsIntoTSVFile)
class AggregateMergedOutput(DefaultPipelineTask):

    def output(self):
        artifacts_dir = self.cfg.output_dir + "/release/artifacts/"
        return luigi.LocalTarget(artifacts_dir + "aggregated.tsv")

    def run(self):
        artifacts_dir = self.cfg.output_dir + "/release/artifacts/"
        os.chdir(data_merging_method_dir)

        args = ["python", "aggregate_across_columns.py", "-i",
                artifacts_dir + "merged.tsv",
                "-o", artifacts_dir + "aggregated.tsv"]
        print "Running aggregate_across_columns.py with the following args: %s" % (
            args)
        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_input_and_output_tsvs_for_same_number_variants(
            artifacts_dir + "merged.tsv",
            artifacts_dir + "aggregated.tsv")


@requires(AggregateMergedOutput)
class BuildAggregatedOutput(DefaultPipelineTask):

    def output(self):
        artifacts_dir = self.cfg.output_dir + "/release/artifacts/"
        return luigi.LocalTarget(artifacts_dir + "built.tsv")

    def run(self):
        release_dir = self.cfg.output_dir + "/release/"
        artifacts_dir = release_dir + "artifacts/"
        brca_resources_dir = self.cfg.resources_dir
        os.chdir(data_merging_method_dir)

        args = ["python", "brca_pseudonym_generator.py",
                artifacts_dir + "aggregated.tsv",
                artifacts_dir + "built.tsv",
                "--log-path", artifacts_dir + "brca-pseudonym-generator.log",
                "--config-file", self.cfg.gene_config_path,
                "--resources", brca_resources_dir]

        print "Running brca_pseudonym_generator.py with the following args: %s" % (
            args)

        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_input_and_output_tsvs_for_same_number_variants(
            artifacts_dir + "aggregated.tsv",
            artifacts_dir + "built.tsv")


@requires(BuildAggregatedOutput)
class AppendMupitStructure(DefaultPipelineTask):

    def output(self):
        artifacts_dir = self.cfg.output_dir + "/release/artifacts/"
        return luigi.LocalTarget(artifacts_dir + "built_with_mupit.tsv")

    def run(self):
        release_dir = self.cfg.output_dir + "/release/"
        artifacts_dir = release_dir + "artifacts/"
        brca_resources_dir = self.cfg.resources_dir
        os.chdir(data_merging_method_dir)

        args = ["python", "getMupitStructure.py", "-i",
                artifacts_dir + "built.tsv", "-o",
                artifacts_dir + "/built_with_mupit.tsv"]
        print "Running getMupitStructure.py with the following args: %s" % (
            args)
        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_input_and_output_tsvs_for_same_number_variants(
            artifacts_dir + "built.tsv",
            artifacts_dir + "built_with_mupit.tsv")


@requires(AppendMupitStructure)
class CalculatePriors(DefaultPipelineTask):
    def output(self):
        artifacts_dir = self.cfg.output_dir + "/release/artifacts/"
        return luigi.LocalTarget(artifacts_dir + "built_with_priors.tsv")

    def run(self):
        artifacts_dir_host = self.cfg.output_dir_host + "/release/artifacts/"
        os.chdir(priors_method_dir)

        args = ['bash', 'calcpriors.sh', self.cfg.priors_references_dir,
                artifacts_dir_host, 'built_with_mupit.tsv',
                'built_with_priors.tsv', self.cfg.priors_docker_image_name]

        print "Running calcpriors.sh with the following args: %s" % (args)
        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_input_and_output_tsvs_for_same_number_variants(
            self.input().path,
            self.output().path)


@requires(CalculatePriors)
class FilterBlacklistedPriors(DefaultPipelineTask):
    def output(self):
        artifacts_dir = self.cfg.output_dir + "/release/artifacts/"
        return luigi.LocalTarget(artifacts_dir + "built_with_priors_clean.tsv")

    def run(self):
        artifacts_dir = self.cfg.output_dir + "/release/artifacts/"
        os.chdir(priors_filter_method_dir)

        args = ["python", "filterBlacklistedVars.py",
                "--output", artifacts_dir + "built_with_priors_clean.tsv",
                "--blacklisted_vars", "blacklisted_vars.txt",
                "filter",
                artifacts_dir + "built_with_priors.tsv"]

        print "Running filterBlacklistedVars.py with the following args: %s" % (
            args)
        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        # we only clear a few columns; we shouldn't be gaining or losing any variants
        pipeline_utils.check_input_and_output_tsvs_for_same_number_variants(
            self.input().path,
            self.output().path)


@requires(FilterBlacklistedPriors)
class AppendVRId(DefaultPipelineTask):
    def output(self):
        artifacts_dir = self.cfg.output_dir + "/release/artifacts/"
        return luigi.LocalTarget(artifacts_dir + "built_with_VR_ids.tsv")

    def run(self):
        artifacts_dir_host = self.cfg.output_dir_host + "/release/artifacts/"
        os.chdir(vr_method_dir)

        args = [
            'bash', 'appendvrids.sh',
            artifacts_dir_host,
            'built_with_priors_clean.tsv',
            'built_with_vr_ids.tsv',
            self.cfg.vr_docker_image_name,
            self.cfg.seq_repo_dir
        ]


        print "Running appendVRIds.py with the following args: %s" % (
            args)
        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        # we shouldn't be gaining or losing any variants
        pipeline_utils.check_input_and_output_tsvs_for_same_number_variants(
            self.input().path,
            self.output().path)


@requires(AppendVRId)
class FindMissingReports(DefaultPipelineTask):
    def output(self):
        artifacts_dir = self.cfg.output_dir + "/release/artifacts/"
        return luigi.LocalTarget(artifacts_dir + "missing_reports.log")

    def run(self):
        release_dir = self.cfg.output_dir + "/release/"
        artifacts_dir = self.cfg.output_dir + "/release/artifacts/"
        os.chdir(data_merging_method_dir)

        args = ["python", "check_for_missing_reports.py", "-b",
                artifacts_dir + "built_with_VR_ids.tsv", "-r",
                artifacts_dir,
                "-a", artifacts_dir, "-v"]
        print "Running check_for_missing_reports.py with the following args: %s" % (
            args)
        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(
            artifacts_dir + "missing_reports.log")


@requires(FindMissingReports)
class RunDiffAndAppendChangeTypesToOutput(DefaultPipelineTask):
    def _extract_release_date(self, version_json):
        with open(version_json, 'r') as f:
            j = json.load(f)
            return datetime.datetime.strptime(j['date'], '%Y-%m-%d')

    def output(self):
        release_dir = self.cfg.output_dir + "/release/"
        diff_dir = pipeline_utils.create_path_if_nonexistent(
            release_dir + "diff/")
        return {'built_with_change_types': luigi.LocalTarget(
            release_dir + "built_with_change_types.tsv"),
                'removed': luigi.LocalTarget(diff_dir + "removed.tsv"),
                'added': luigi.LocalTarget(diff_dir + "added.tsv"),
                'added_data': luigi.LocalTarget(diff_dir + "added_data.tsv"),
                'diff': luigi.LocalTarget(diff_dir + "diff.txt"),
                'diff_json': luigi.LocalTarget(diff_dir + "diff.json"),
                'README': luigi.LocalTarget(diff_dir + "README.txt")}

    def run(self):
        release_dir = self.cfg.output_dir + "/release/"
        artifacts_dir = release_dir + "artifacts/"
        diff_dir = pipeline_utils.create_path_if_nonexistent(
            release_dir + "diff/")
        os.chdir(utilities_method_dir)

        tmp_dir = tempfile.mkdtemp()
        previous_data_path = pipeline_utils.extract_file(
            self.cfg.previous_release_tar, tmp_dir,
            'output/release/built_with_change_types.tsv')
        version_json_path = pipeline_utils.extract_file(
            self.cfg.previous_release_tar, tmp_dir,
            'output/release/metadata/version.json')
        previous_release_date = self._extract_release_date(version_json_path)
        previous_release_date_str = datetime.datetime.strftime(
            previous_release_date, '%m-%d-%Y')

        args = ["python", "releaseDiff.py", "--v2",
                artifacts_dir + "built_with_VR_ids.tsv", "--v1",
                previous_data_path,
                "--removed", diff_dir + "removed.tsv", "--added",
                diff_dir + "added.tsv", "--added_data",
                diff_dir + "added_data.tsv", "--diff", diff_dir + "diff.txt",
                "--diff_json", diff_dir + "diff.json",
                "--output", release_dir + "built_with_change_types.tsv",
                "--artifacts_dir", artifacts_dir,
                "--diff_dir", diff_dir, "--v1_release_date",
                previous_release_date_str, "--reports", "False"]

        print "Running releaseDiff.py with the following args: %s" % (args)
        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        shutil.rmtree(tmp_dir)  # cleaning up

        pipeline_utils.check_input_and_output_tsvs_for_same_number_variants(
            artifacts_dir + "built_with_VR_ids.tsv",
            release_dir + "built_with_change_types.tsv")


@requires(RunDiffAndAppendChangeTypesToOutput)
class RunDiffAndAppendChangeTypesToOutputReports(DefaultPipelineTask):
    def _extract_release_date(self, version_json):
        with open(version_json, 'r') as f:
            j = json.load(f)
            return datetime.datetime.strptime(j['date'], '%Y-%m-%d')

    def output(self):
        release_dir = self.cfg.output_dir + "/release/"
        diff_dir = pipeline_utils.create_path_if_nonexistent(
            release_dir + "diff/")
        return {'reports_with_change_types': luigi.LocalTarget(
            release_dir + "reports_with_change_types.tsv"),
                'removed_reports': luigi.LocalTarget(
                    diff_dir + "removed_reports.tsv"),
                'added_reports': luigi.LocalTarget(
                    diff_dir + "added_reports.tsv"),
                'added_data_reports': luigi.LocalTarget(
                    diff_dir + "added_data_reports.tsv"),
                'diff_reports': luigi.LocalTarget(
                    diff_dir + "diff_reports.txt"),
                'diff_json_reports': luigi.LocalTarget(
                    diff_dir + "diff_reports.json"),
                'README': luigi.LocalTarget(diff_dir + "README.txt")}

    def run(self):
        release_dir = self.cfg.output_dir + "/release/"
        artifacts_dir = release_dir + "artifacts/"
        diff_dir = pipeline_utils.create_path_if_nonexistent(
            release_dir + "diff/")
        os.chdir(utilities_method_dir)

        tmp_dir = tempfile.mkdtemp()
        previous_data_path = pipeline_utils.extract_file(
            self.cfg.previous_release_tar, tmp_dir,
            'output/release/artifacts/reports.tsv')
        version_json_path = pipeline_utils.extract_file(
            self.cfg.previous_release_tar, tmp_dir,
            'output/release/metadata/version.json')
        previous_release_date = self._extract_release_date(version_json_path)
        previous_release_date_str = datetime.datetime.strftime(
            previous_release_date, '%m-%d-%Y')

        args = ["python", "releaseDiff.py", "--v2",
                artifacts_dir + "reports.tsv", "--v1", previous_data_path,
                "--removed", diff_dir + "removed_reports.tsv", "--added",
                diff_dir + "added_reports.tsv", "--added_data",
                diff_dir + "added_data_reports.tsv", "--diff",
                diff_dir + "diff_reports.txt", "--diff_json",
                diff_dir + "diff_reports.json",
                "--output", release_dir + "reports_with_change_types.tsv",
                "--artifacts_dir", artifacts_dir,
                "--diff_dir", diff_dir, "--v1_release_date",
                previous_release_date_str, "--reports", "True"]

        print "Running releaseDiff.py with the following args: %s" % (args)
        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        shutil.rmtree(tmp_dir)  # cleaning up

        pipeline_utils.check_input_and_output_tsvs_for_same_number_variants(
            artifacts_dir + "reports.tsv",
            release_dir + "reports_with_change_types.tsv")


@requires(RunDiffAndAppendChangeTypesToOutputReports)
class GenerateReleaseNotes(DefaultPipelineTask):

    def output(self):
        metadata_dir = pipeline_utils.create_path_if_nonexistent(
            self.cfg.output_dir + "/release/metadata/")
        return luigi.LocalTarget(metadata_dir + "version.json")

    def run(self):
        metadata_dir = self.cfg.output_dir + "/release/metadata/"
        os.chdir(data_merging_method_dir)

        args = ["python", "buildVersionMetadata.py", "--date",
                str(self.cfg.date), "--notes", self.cfg.release_notes,
                "--output", metadata_dir + "version.json"]
        print "Running buildVersionMetadata.py with the following args: %s" % (
            args)
        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(metadata_dir + "version.json")


@requires(GenerateReleaseNotes)
class TopLevelReadme(DefaultPipelineTask):
    def output(self):
        top_level_readme_dest = os.path.join(self.cfg.output_dir, "README.txt")
        return luigi.LocalTarget(top_level_readme_dest)

    def run(self):
        top_level_readme_src = os.path.abspath(
            os.path.join(os.path.realpath(__file__), os.pardir, os.pardir,
                         "top_level_readme.txt"))

        shutil.copyfile(top_level_readme_src, self.output().path)


@requires(TopLevelReadme)
class DataDictionary(DefaultPipelineTask):
    def output(self):
        release_dir = self.cfg.output_dir + "/release/"
        data_dictionary_dest = os.path.join(release_dir, "built_with_change_types.dictionary.tsv")
        return luigi.LocalTarget(data_dictionary_dest)

    def run(self):
        data_dictionary_src = os.path.abspath(
            os.path.join(os.path.realpath(__file__), os.pardir, os.pardir, "built_with_change_types.dictionary.tsv"))

        shutil.copyfile(data_dictionary_src, self.output().path)


@requires(DataDictionary)
class GenerateMD5Sums(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(self.cfg.output_dir + "/md5sums.txt")

    def run(self):
        output_dir = self.cfg.output_dir
        md5sumsFile = output_dir + "/md5sums.txt"

        os.chdir(utilities_method_dir)

        args = ["python", "generateMD5Sums.py", "-i", output_dir, "-o",
                md5sumsFile]
        print "Generating md5sums with the following args: %s" % (args)
        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(md5sumsFile)


@requires(GenerateMD5Sums)
class GenerateReleaseArchive(DefaultPipelineTask):

    def getArchiveName(self):
        # Format archive filename as release-mm-dd-yy.tar.gz
        return "release-" + self.cfg.date.strftime("%x").replace('/',
                                                                 '-') + ".tar.gz"

    def getArchiveParentDirectory(self):
        return os.path.dirname(self.cfg.output_dir) + "/"

    def output(self):
        return luigi.LocalTarget(
            self.getArchiveParentDirectory() + self.getArchiveName())

    def run(self):
        os.chdir(self.getArchiveParentDirectory())
        with tarfile.open(
                self.getArchiveParentDirectory() + self.getArchiveName(),
                "w:gz") as tar:
            tar.add(self.cfg.output_dir,
                    arcname=os.path.basename(self.cfg.output_dir))
