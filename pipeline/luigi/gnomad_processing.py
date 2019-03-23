import os
import subprocess
import logging

import luigi
from luigi.util import requires

import pipeline_utils
from pipeline_common import DefaultPipelineTask

###############################################
#                  gnomAD                     #
###############################################

gnomAD_method_dir = os.path.abspath('../gnomad')
gnomAD_brca1_data_url = "https://brcaexchange.org/backend/downloads/BRCA1_gnomAD_v2.1.1_(non-cancer)_ENSG00000012048_2019_03_21_12_39_22.csv"
gnomAD_brca2_data_url = "https://brcaexchange.org/backend/downloads/BRCA2_gnomAD_v2.1.1_(non-cancer)_ENSG00000139618_2019_03_21_12_40_20.csv"

logger = logging.getLogger('gnomAD')


class GnomADTask(DefaultPipelineTask):
    def __init__(self, *args, **kwargs):
        super(GnomADTask, self).__init__(*args, **kwargs)
        self.gnomAD_file_dir = self.cfg.file_parent_dir + "/gnomAD"
        self.gnomAD_brca1_file_name = gnomAD_brca1_data_url.split('/')[-1]
        self.gnomAD_brca2_file_name = gnomAD_brca2_data_url.split('/')[-1]

        pipeline_utils.create_path_if_nonexistent(self.gnomAD_file_dir)


class DownloadBRCA1GnomADData(GnomADTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.gnomAD_file_dir, self.gnomAD_brca1_file_name))

    def run(self):
        os.chdir(self.gnomAD_file_dir)

        pipeline_utils.download_file_and_display_progress(gnomAD_brca1_data_url)


@requires(DownloadBRCA1GnomADData)
class DownloadBRCA2GnomADData(GnomADTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.gnomAD_file_dir, self.gnomAD_brca2_file_name))

    def run(self):
        os.chdir(self.gnomAD_file_dir)

        pipeline_utils.download_file_and_display_progress(gnomAD_brca2_data_url)


@requires(DownloadBRCA2GnomADData)
class ParseGnomADBRCA1Data(GnomADTask):

    def output(self):
        return luigi.LocalTarget(os.path.join(self.gnomAD_file_dir, "gnomAD_BRCA1.clean.tsv"))

    def run(self):
        os.chdir(gnomAD_method_dir)

        args = ["python", "parse_gnomad_data.py", "-i", os.path.join(self.gnomAD_file_dir, self.gnomAD_brca1_file_name), "-o",
                self.output().path]

        logger.info("Running parse_gnomad_data.py with the following args: %s", args)

        sp = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(self.output().path)


@requires(ParseGnomADBRCA1Data)
class ParseGnomADBRCA2Data(GnomADTask):

    def output(self):
        return luigi.LocalTarget(os.path.join(self.gnomAD_file_dir, "gnomAD_BRCA2.clean.tsv"))

    def run(self):
        os.chdir(gnomAD_method_dir)

        args = ["python", "parse_gnomad_data.py", "-i", os.path.join(self.gnomAD_file_dir, self.gnomAD_brca2_file_name), "-o",
                self.output().path]

        logger.info("Running parse_gnomad_data.py with the following args: %s", args)

        sp = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(self.output().path)



@requires(ParseGnomADBRCA2Data)
class ConvertGnomADBRCA1ToVCF(GnomADTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.gnomAD_file_dir, 'gnomAD_BRCA1.clean.hg19.vcf'))

    def run(self):
        brca_resources_dir = self.cfg.resources_dir
        artifacts_dir = pipeline_utils.create_path_if_nonexistent(self.cfg.output_dir + "/release/artifacts")

        os.chdir(gnomAD_method_dir)

        args = ["python", "gnomad_to_vcf.py", "-i", os.path.join(self.gnomAD_file_dir, "gnomAD_BRCA1.clean.tsv"), "-o",
                self.output().path, "-a", "gnomADAnnotation",
                "-r", brca_resources_dir + "/refseq_annotation.hg19.gp", "-g",
                brca_resources_dir + "/hg19.fa", "-l", artifacts_dir + "/gnomAD_BRCA1_error_variants.log",
                "-s", "gnomAD"]

        logger.info("Running gnomad_to_vcf.py with the following args: %s", args)

        sp = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(self.output().path)


@requires(ConvertGnomADBRCA1ToVCF)
class ConvertGnomADBRCA2ToVCF(GnomADTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.gnomAD_file_dir, 'gnomAD_BRCA2.clean.hg19.vcf'))

    def run(self):
        brca_resources_dir = self.cfg.resources_dir
        artifacts_dir = pipeline_utils.create_path_if_nonexistent(self.cfg.output_dir + "/release/artifacts")

        os.chdir(gnomAD_method_dir)

        args = ["python", "gnomad_to_vcf.py", "-i", os.path.join(self.gnomAD_file_dir, "gnomAD_BRCA2.clean.tsv"), "-o",
                self.output().path, "-a", "gnomADAnnotation",
                "-r", brca_resources_dir + "/refseq_annotation.hg19.gp", "-g",
                brca_resources_dir + "/hg19.fa", "-l", artifacts_dir + "/gnomAD_BRCA2_error_variants.log",
                "-s", "gnomAD"]

        logger.info("Running gnomad_to_vcf.py with the following args: %s", args)

        sp = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(self.output().path)


@requires(ConvertGnomADBRCA2ToVCF)
class ConcatenateBRCA12GnomADData(GnomADTask):

    def output(self):
        return luigi.LocalTarget(os.path.join(self.gnomAD_file_dir, "gnomAD_BRCA12.clean.hg19.vcf"))

    def run(self):
        concatenated_brca_output_file = self.output().path

        with open(concatenated_brca_output_file, 'w') as f:

            args = ["vcf-concat", self.gnomAD_file_dir + "/gnomAD_BRCA1.clean.hg19.vcf",
                    self.gnomAD_file_dir + "/gnomAD_BRCA2.clean.hg19.vcf"]

            logger.info("Calling vcf-concat with the following args: %s", args)

            sp = subprocess.Popen(args, stdout=f, stderr=subprocess.PIPE)
            pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(self.output().path)
        logger.info("Concatenation complete.")



@requires(ConcatenateBRCA12GnomADData)
class CrossmapGnomADBRCA12Data(GnomADTask):

    def output(self):
        return luigi.LocalTarget(os.path.join(self.gnomAD_file_dir, "gnomAD_BRCA12.clean.hg38.vcf"))

    def run(self):
        brca_resources_dir = self.cfg.resources_dir

        args = ["CrossMap.py", "vcf", brca_resources_dir + "/hg19ToHg38.over.chain.gz",
                self.gnomAD_file_dir + "/gnomAD_BRCA12.clean.hg19.vcf", brca_resources_dir + "/hg38.fa",
                self.gnomAD_file_dir + "/gnomAD_BRCA12.clean.hg38.vcf"]

        logger.info("Running CrossMap.py with the following args: %s", args)

        sp = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(self.gnomAD_file_dir + "/gnomAD_BRCA12.clean.hg38.vcf")


@requires(CrossmapGnomADBRCA12Data)
class SortConcatenatedGnomADData(GnomADTask):

    def output(self):
        return luigi.LocalTarget(os.path.join(self.gnomAD_file_dir, "gnomad_BRCA12.clean.sorted.hg38.vcf"))

    def run(self):
        concatenated_brca_output_file = self.input().path
        sorted_concatenated_brca_output_file = self.output().path

        with open(concatenated_brca_output_file, 'w') as f:
            args = ["vcf-sort", self.gnomAD_file_dir + "/gnomad_BRCA12.clean.hg38.vcf"]
            logger.info("Running vcf-sort with the following args: %s", args)
            sp = subprocess.Popen(args, stdout=f, stderr=subprocess.PIPE)
            pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(sorted_concatenated_brca_output_file)
        logger.info("Sorting of concatenated files complete.")
