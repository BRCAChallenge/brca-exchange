import logging
import os
import luigi
from luigi.util import requires

from workflow import pipeline_utils
from workflow.pipeline_common import DefaultPipelineTask

import pdb

###############################################
#                  gnomAD                     #
###############################################

gnomAD_method_dir = os.path.abspath('../gnomad')

logger = logging.getLogger('gnomAD')


class GnomADTask(DefaultPipelineTask):
    def __init__(self, *args, **kwargs):
        super(GnomADTask, self).__init__(*args, **kwargs)

    def _download_file(self, url, output_path):
        data = pipeline_utils.urlopen_with_retry(url).read()
        with open(output_path, "wb") as f:
            f.write(data)


class DownloadGnomADData(GnomADTask):
    gnomAD_v2_static_data_url = luigi.Parameter(default='https://brcaexchange.org/backend/downloads/gnomAD_v2_hg19_10_02_2020.tsv',
                                            description='URL to download static gnomAD v2 data from')
    gnomAD_v3_static_data_url = luigi.Parameter(default='https://brcaexchange.org/backend/downloads/gnomAD_v3_1_1_GRCh38_06_07_2021.tsv',
                                            description='URL to download static gnomAD v3 data from')

    def output(self):
        return { "v2": luigi.LocalTarget(f"{self.gnomad_file_dir}/gnomAD_v2_hg19_10_02_2020.tsv"),
                 "v3": luigi.LocalTarget(f"{self.gnomad_file_dir}/gnomAD_v3_1_1_GRCh38_06_07_2021.tsv")}


    def run(self):
        self._download_file(self.gnomAD_v2_static_data_url, self.output()["v2"].path)
        self._download_file(self.gnomAD_v3_static_data_url, self.output()["v3"].path)


@requires(DownloadGnomADData)
class ConvertGnomADToVCF(GnomADTask):
    def output(self):
        return { "v2": luigi.LocalTarget(f"{self.gnomad_file_dir}/gnomADv2.hg19.vcf"),
                 "v3": luigi.LocalTarget(f"{self.gnomad_file_dir}/gnomADv3.hg38.vcf")}

    def run(self):
        os.chdir(gnomAD_method_dir)

        for k, v in self.input().items():
            if k == "v3":
                annotation = "gnomADv3Annotation"
            else:
                annotation = "gnomADAnnotation"
            args = ["python", "gnomad_to_vcf.py", "-i", v.path, "-o",
                    self.output()[k].path, "-a", annotation,
                    "-l", self.artifacts_dir + "/gnomADv2_error_variants.log",
                    "-s", "gnomAD"]

            pipeline_utils.run_process(args)
            pipeline_utils.check_file_for_contents(self.output()[k].path)


@requires(ConvertGnomADToVCF)
class CrossmapGnomADV2Data(GnomADTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.gnomad_file_dir, "gnomADv2.hg38.vcf"))

    def run(self):
        brca_resources_dir = self.cfg.resources_dir

        args = ["CrossMap.py", "vcf", brca_resources_dir + "/hg19ToHg38.over.chain.gz",
                self.input()["v2"].path, brca_resources_dir + "/hg38.fa",
                self.output().path]

        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(CrossmapGnomADV2Data)
class SortGnomADData(GnomADTask):
    def output(self):
        return { "v2": luigi.LocalTarget(f"{self.gnomad_file_dir}/gnomADv2.sorted.hg38.vcf"),
                 "v3": luigi.LocalTarget(f"{self.gnomad_file_dir}/gnomADv3.sorted.hg38.vcf")}

    def run(self):
        for key in self.output().keys():
            args = ["vcf-sort", f"{self.gnomad_file_dir}/gnomAD{key}.hg38.vcf"]
            pipeline_utils.run_process(args, redirect_stdout_path=self.output()[key].path)
            pipeline_utils.check_file_for_contents(self.output()[key].path)

            
class DownloadStaticGnomADVCF(GnomADTask):
    gnomAD_v2_static_vcf_url = luigi.Parameter(default='https://brcaexchange.org/backend/downloads/gnomADv2.sorted.hg38.vcf',
                                               description='URL to download static gnomAD v2 vcf from')
    gnomAD_v3_static_data_url = luigi.Parameter(default='https://brcaexchange.org/backend/downloads/gnomADv3.sorted.hg38.vcf',
                                            description='URL to download static gnomAD v3 vcf from')
    
    def output(self):
        return { "v2": luigi.LocalTarget(f"{self.gnomad_file_dir}/gnomADv2.sorted.hg38.vcf"),
                 "v3": luigi.LocalTarget(f"{self.gnomad_file_dir}/gnomADv3.sorted.hg38.vcf")}
    

    def run(self):
        self._download_file(self.gnomAD_v2_static_vcf_url, self.output()["v2"].path)
        self._download_file(self.gnomAD_v3_static_data_url, self.output()["v3"].path)


class DownloadGnomADv4(GnomADTask):
    def output(self):
        return luigi.LocalTarget(f"{self.gnomad_file_dir}/gnomADv4.sorted.hg38.vcf") 

    def run(self):
        os.chdir(gnomAD_method_dir)
        args = ["python", "download_gnomad_fourpointone.py",
                "-g", self.cfg.gene_config_path,
                "-o", self.output().path]
        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)
