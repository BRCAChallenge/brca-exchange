import os
import subprocess
import tarfile
import logging

import luigi
from luigi.util import requires

import pipeline_utils
from pipeline_common import DefaultPipelineTask

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

        logger.info(
            "Calling espExtract.py for %s region with the following arguments: %s",
            self.gene,
            args)
        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)


class ConcatenateESPData(ESPTask):
    def requires(self):
        for gene in self.cfg.gene_metadata['symbol']:
            yield ExtractESPData(gene)

    def output(self):
        return luigi.LocalTarget(self.esp_file_dir + "/esp.hg38.vcf")

    def run(self):
        concatenated_brca_output_file = self.output().path

        with open(concatenated_brca_output_file, 'w') as f:
            # Note: requires correct installation of VCF tools and export PERL5LIB=/path/to/your/vcftools-directory/src/perl/ in path
            args = ["vcf-concat"]
            args.extend([t.path for t in self.input()])

            logger.info("Calling vcf-concat with the following args: %s", args)
            sp = subprocess.Popen(args, stdout=f, stderr=subprocess.PIPE)
            pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(concatenated_brca_output_file)
        logger.info("Concatenation complete.")


@requires(ConcatenateESPData)
class SortConcatenatedESPData(ESPTask):
    def output(self):
        return luigi.LocalTarget(self.esp_file_dir + "/esp.sorted.hg38.vcf")

    def run(self):
        sorted_concatenated_brca_output_file = self.output().path
        concatenated_brca_output_file = self.input().path

        with open(sorted_concatenated_brca_output_file, 'w') as f:
            args = ["vcf-sort", concatenated_brca_output_file]
            logger.info("Calling vcf-sort with the following args: %s", args)
            sp = subprocess.Popen(args, stdout=f, stderr=subprocess.PIPE)
            pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_file_for_contents(
            sorted_concatenated_brca_output_file)
        logger.info("Sorting of concatenated files complete.")
