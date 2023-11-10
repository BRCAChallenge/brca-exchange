import logging
import os
import tarfile

import luigi
from luigi.util import requires

from workflow import pipeline_utils
from workflow.pipeline_common import DefaultPipelineTask

###############################################
#                     ESP                     #
###############################################

esp_method_dir = os.path.abspath('../esp')
esp_data_url = "http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.tar.gz"

logger = logging.getLogger('ESP')


class ESPTask(DefaultPipelineTask):
    def __init__(self, *args, **kwargs):
        super(ESPTask, self).__init__(*args, **kwargs)
        self.esp_file_dir = self.cfg.file_parent_dir + "/ESP"
        self.esp_tar_file_name = esp_data_url.split('/')[-1]

        pipeline_utils.create_path_if_nonexistent(self.esp_file_dir)


class DownloadStaticESPData(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(f"{self.esp_file_dir}/esp.sorted.hg38.vcf")

    def run(self):
        os.chdir(self.esp_file_dir)

        esp_vcf_url = "https://brcaexchange.org/backend/downloads/esp.sorted.hg38.vcf"
        pipeline_utils.download_file_and_display_progress(esp_vcf_url)


        
class DownloadLatestESPData(ESPTask):
    def output(self):
        return luigi.LocalTarget(
            os.path.join(self.esp_file_dir, self.esp_tar_file_name))

    def run(self):
        os.chdir(self.esp_file_dir)

        pipeline_utils.download_file_and_display_progress(esp_data_url)


@requires(DownloadLatestESPData)
class DecompressESPTarfile(ESPTask):
    chr = luigi.IntParameter()

    def __init__(self, chr):
        super(DecompressESPTarfile, self).__init__(chr=chr)

    def output(self):
        return luigi.LocalTarget(
            self.esp_file_dir + "/ESP6500SI-V2-SSA137.GRCh38-liftover.chr" + str(
                self.chr) + ".snps_indels.vcf")

    def run(self):
        os.chdir(self.esp_file_dir)

        tar = tarfile.open(self.esp_tar_file_name, "r:gz")
        tar.extract(os.path.basename(self.output().path))
        tar.close()
        logger.info("Finished extracting files from %s", self.esp_tar_file_name)


class ExtractESPData(ESPTask):
    gene = luigi.Parameter()

    def __init__(self, gene):
        super(ExtractESPData, self).__init__(gene)

    def requires(self):
        chr = self.cfg.gene_metadata.loc[self.gene]['chr']
        yield DecompressESPTarfile(chr=chr)

    def output(self):
        return luigi.LocalTarget(
            self.esp_file_dir + "/esp." + self.gene.lower() + ".vcf")

    def run(self):
        gene_record = self.cfg.gene_metadata.loc[self.gene]

        chr_vcf_path = self.input()[0].path

        args = ["python", os.path.join(esp_method_dir, "espExtract.py"),
                chr_vcf_path,
                "--start", str(gene_record['start_hg38']),
                "--end", str(gene_record['end_hg38']),
                "--full", "1", "-o",
                self.output().path]

        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)


class ConcatenateESPData(ESPTask):
    def requires(self):
        for gene in self.cfg.gene_metadata['symbol']:
            yield ExtractESPData(gene)

    def output(self):
        return luigi.LocalTarget(self.esp_file_dir + "/esp.hg38.vcf")

    def run(self):
        # Note: requires correct installation of VCF tools and export PERL5LIB=/path/to/your/vcftools-directory/src/perl/ in path
        args = ["vcf-concat"]
        args.extend([t.path for t in self.input()])

        pipeline_utils.run_process(args, redirect_stdout_path=self.output().path)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(ConcatenateESPData)
class DownloadSortedESPData(ESPTask):
    def output(self):
        return luigi.LocalTarget(self.esp_file_dir + "/esp.sorted.hg38.vcf")

    def run(self):
        args = ["vcf-sort", self.input().path]

        pipeline_utils.run_process(args, redirect_stdout_path=self.output().path)
        pipeline_utils.check_file_for_contents(self.output().path)
