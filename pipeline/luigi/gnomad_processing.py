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

logger = logging.getLogger('gnomAD')


class GnomADTask(DefaultPipelineTask):
    def __init__(self, *args, **kwargs):
        super(GnomADTask, self).__init__(*args, **kwargs)
        self.gnomAD_file_dir = self.cfg.file_parent_dir + "/gnomAD"
        self.gnomAD_static_file = "gnomAD.concatted.clean.tsv"

        pipeline_utils.create_path_if_nonexistent(self.gnomAD_file_dir)


"""
Until issues are resolved with the GnomAD API, use static concattenated data

class DownloadGnomADData(GnomADTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.gnomAD_file_dir, self.gnomAD_download_file_name))

    def run(self):
        artifacts_dir = pipeline_utils.create_path_if_nonexistent(self.cfg.output_dir + "/release/artifacts")
        os.chdir(gnomAD_method_dir)

        args = ["python", "download_gnomad_data.py", "-o", self.output().path,
                "-l", artifacts_dir + "/download_gnomAD_data.log"]
        sp = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(self.output().path)


@requires(DownloadGnomADData)
class NormalizeGnomADData(GnomADTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.gnomAD_file_dir, "gnomAD.normalized.tsv"))

    def run(self):
        os.chdir(gnomAD_method_dir)

        args = ["python", "normalize.py", "-i", self.input().path,
                "-o", self.output().path]

        logger.info("Running normalize.py with the following args: %s", args)

        sp = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(self.output().path)


@requires(NormalizeGnomADData)
class ParseGnomADData(GnomADTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.gnomAD_file_dir, "gnomAD.clean.tsv"))

    def run(self):
        os.chdir(gnomAD_method_dir)

        args = ["python", "parse_gnomad_data.py", "-i", self.input().path,
                "-o", self.output().path]

        logger.info("Running parse_gnomad_data.py with the following args: %s", args)

        sp = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(self.output().path)
"""

class DownloadStaticGnomADData(GnomADTask):
    def output(self):
        return luigi.LocalTarget(
            os.path.join(self.gnomAD_file_dir, self.gnomAD_static_file))

    def run(self):
        os.chdir(self.gnomAD_file_dir)
        static_gnomAD_data_url = "https://brcaexchange.org/backend/downloads/gnomAD.concatted.clean.tsv"
        pipeline_utils.download_file_and_display_progress(static_gnomAD_data_url, self.output().path)


@requires(DownloadStaticGnomADData)
class ConvertGnomADToVCF(GnomADTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.gnomAD_file_dir, 'gnomAD.clean.hg19.vcf'))

    def run(self):
        brca_resources_dir = self.cfg.resources_dir
        artifacts_dir = pipeline_utils.create_path_if_nonexistent(self.cfg.output_dir + "/release/artifacts")

        os.chdir(gnomAD_method_dir)

        args = ["python", "gnomad_to_vcf.py", "-i", self.input().path, "-o",
                self.output().path, "-a", "gnomADAnnotation",
                "-r", brca_resources_dir + "/refseq_annotation.hg19.gp", "-g",
                brca_resources_dir + "/hg19.fa", "-l", artifacts_dir + "/gnomAD_error_variants.log",
                "-s", "gnomAD"]

        logger.info("Running gnomad_to_vcf.py with the following args: %s", args)

        sp = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(self.output().path)


@requires(ConvertGnomADToVCF)
class CrossmapGnomADData(GnomADTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.gnomAD_file_dir, "gnomAD.clean.hg38.vcf"))

    def run(self):
        brca_resources_dir = self.cfg.resources_dir

        args = ["CrossMap.py", "vcf", brca_resources_dir + "/hg19ToHg38.over.chain.gz",
                self.gnomAD_file_dir + "/gnomAD.clean.hg19.vcf", brca_resources_dir + "/hg38.fa",
                self.gnomAD_file_dir + "/gnomAD.clean.hg38.vcf"]

        logger.info("Running CrossMap.py with the following args: %s", args)

        sp = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(self.output().path)


@requires(CrossmapGnomADData)
class SortGnomADData(GnomADTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.gnomAD_file_dir, "gnomAD.clean.sorted.hg38.vcf"))

    def run(self):
        with open(self.output().path, 'w') as f:
            args = ["vcf-sort", self.input().path]
            logger.info("Running vcf-sort with the following args: %s", args)
            sp = subprocess.Popen(args, stdout=f, stderr=subprocess.PIPE)
            pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(self.output().path)
        logger.info("Sorting of concatenated files complete.")
