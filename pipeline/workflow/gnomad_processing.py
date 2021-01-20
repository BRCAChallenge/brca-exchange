import logging
import os

import luigi
from luigi.util import requires

from workflow import pipeline_utils
from workflow.pipeline_common import DefaultPipelineTask

###############################################
#                  gnomAD                     #
###############################################

gnomAD_method_dir = os.path.abspath('../gnomad')

logger = logging.getLogger('gnomAD')


class GnomADTask(DefaultPipelineTask):
    def __init__(self, *args, **kwargs):
        super(GnomADTask, self).__init__(*args, **kwargs)
        self.gnomAD_download_file = "gnomAD.tsv"


class DownloadGnomADData(GnomADTask):
    gnomAD_static_data_url = luigi.Parameter(default='https://brcaexchange.org/backend/downloads/gnomAD_static_download_10_02_2020.tsv',
                                            description='URL to download static gnomAD data from')

    def output(self):
        return luigi.LocalTarget(self.assays_dir + "/gnomAD.tsv")

    def run(self):
        data = pipeline_utils.urlopen_with_retry(self.gnomAD_static_data_url).read()
        with open(self.output().path, "wb") as f:
            f.write(data)

"""
################
NOTE:
gnomAD rarely updates its dataset
due to issues with consistency downloading the same gnomad data for each release,
a static file can be reused until new data is available
see previous task for static data download
################


class DownloadGnomADData(GnomADTask):
    def output(self):
        return luigi.LocalTarget(
            os.path.join(self.gnomad_file_dir, self.gnomAD_download_file))

    def run(self):
        os.chdir(gnomAD_method_dir)

        args = ["python", "download_gnomad_data.py", "-o", self.output().path, "-l", self.artifacts_dir + "/gnomAD_download.log"]

        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)
"""


@requires(DownloadGnomADData)
class ConvertGnomADToVCF(GnomADTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.gnomad_file_dir, 'gnomAD.hg19.vcf'))

    def run(self):
        os.chdir(gnomAD_method_dir)

        args = ["python", "gnomad_to_vcf.py", "-i", self.input().path, "-o",
                self.output().path, "-a", "gnomADAnnotation",
                "-l", self.artifacts_dir + "/gnomAD_error_variants.log",
                "-s", "gnomAD"]

        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(ConvertGnomADToVCF)
class CrossmapGnomADData(GnomADTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.gnomad_file_dir, "gnomAD.hg38.vcf"))

    def run(self):
        brca_resources_dir = self.cfg.resources_dir

        args = ["CrossMap.py", "vcf", brca_resources_dir + "/hg19ToHg38.over.chain.gz",
                self.input().path, brca_resources_dir + "/hg38.fa",
                self.output().path]

        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(CrossmapGnomADData)
class SortGnomADData(GnomADTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.gnomad_file_dir, "gnomAD.sorted.hg38.vcf"))

    def run(self):
        args = ["vcf-sort", self.input().path]
        pipeline_utils.run_process(args, redirect_stdout_path=self.output().path)
        pipeline_utils.check_file_for_contents(self.output().path)
