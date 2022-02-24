import datetime
import json
import os
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path
from shutil import copy

import luigi
from luigi.util import requires

luigi.auto_namespace(scope=__name__)

from workflow import bayesdel_processing, esp_processing, gnomad_processing, pipeline_common, pipeline_utils
from workflow.pipeline_common import DefaultPipelineTask, clinvar_method_dir, lovd_method_dir, \
    functional_assays_method_dir, data_merging_method_dir, priors_method_dir, priors_filter_method_dir, \
    utilities_method_dir, vr_method_dir, field_metadata_path, field_metadata_path_additional

from common import utils
from data_merging import generate_variants_output_file

#######################################
# Default Globals / Env / Directories #
#######################################

luigi_dir = os.getcwd()

###############################################
#                   CLINVAR                   #
###############################################

class DownloadLatestClinvarData(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(self.clinvar_file_dir + "/ClinVarFullRelease_00-latest.xml.gz")

    def run(self):
        os.chdir(self.clinvar_file_dir)

        clinvar_data_url = "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz"
        pipeline_utils.download_file_and_display_progress(clinvar_data_url)


@requires(DownloadLatestClinvarData)
class ConvertLatestClinvarDataToXML(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(self.clinvar_file_dir + "/ClinVar.xml")

    def run(self):
        os.chdir(clinvar_method_dir)
        genes_opts = [ s for g in self.cfg.gene_metadata['symbol'] for s in ['--gene', g]]

        pipeline_utils.run_process(["python", "filter_clinvar.py", self.input().path,
                                    self.output().path] + genes_opts)

        pipeline_utils.check_file_for_contents(self.output().path)


@requires(ConvertLatestClinvarDataToXML)
class ConvertClinvarXMLToTXT(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(self.clinvar_file_dir + "/ClinVar.txt")

    def run(self):
        os.chdir(clinvar_method_dir)

        args = ["python", "clinVarParse.py",
                self.input().path,
                "--logs", self.clinvar_file_dir + "/clinvar_xml_to_txt.log",
                "--assembly", "GRCh38"]

        pipeline_utils.run_process(args, redirect_stdout_path=self.output().path)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(ConvertClinvarXMLToTXT)
class ConvertClinvarTXTToVCF(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(self.clinvar_file_dir + "/ClinVar.vcf")

    def run(self):
        os.chdir(data_merging_method_dir)
        args = ["python", "convert_tsv_to_vcf.py", "-i",
                self.input().path, "-o",
                self.output().path, "-s", "ClinVar"]
        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(ConvertClinvarTXTToVCF)
class CopyClinvarVCFToOutputDir(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(self.cfg.output_dir + "/ClinVar.vcf")

    def run(self):
        pipeline_utils.create_path_if_nonexistent(self.cfg.output_dir)

        copy(self.input().path, self.cfg.output_dir)
        pipeline_utils.check_file_for_contents(self.output().path)


###############################################
#                     BIC                     #
###############################################


class DownloadBICData(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(self.bic_file_dir + "/bic_brca12.sorted.hg38.vcf")

    def run(self):
        os.chdir(self.bic_file_dir)

        brca1_data_url = "https://brcaexchange.org/backend/downloads/bic_brca12.sorted.hg38.vcf"
        pipeline_utils.download_file_and_display_progress(brca1_data_url)


@requires(DownloadBICData)
class CopyBICOutputToOutputDir(DefaultPipelineTask):

    def output(self):
        return luigi.LocalTarget(self.cfg.output_dir + "/bic_brca12.sorted.hg38.vcf")

    def run(self):
        copy(self.bic_file_dir + "/bic_brca12.sorted.hg38.vcf", self.cfg.output_dir)
        pipeline_utils.check_file_for_contents(self.output().path)


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
        # calculating host path because we are running a docker within a docker
        ex_lovd_file_dir_host = os.path.join(os.path.dirname(self.cfg.output_dir_host), ExtractDataFromLatestEXLOVD.dir_name)

        os.chdir(lovd_method_dir)

        ex_lovd_data_host_url = "http://hci-exlovd.hci.utah.edu/"

        args = ['bash', 'extract_latest_exlovd.sh', ex_lovd_file_dir_host, "-u", ex_lovd_data_host_url, "-l", "BRCA1",
                "BRCA2", "-o", "/data"]

        pipeline_utils.run_process(args)


@requires(ExtractDataFromLatestEXLOVD)
class ConvertEXLOVDBRCA1ExtractToVCF(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.ex_lovd_file_dir, "exLOVD_brca1.hg19.vcf"))

    def run(self):
        os.chdir(lovd_method_dir)

        args = ["./lovd2vcf.py", "-i", self.ex_lovd_file_dir + "/BRCA1.txt", "-o",
                self.output().path, "-a",
                "exLOVDAnnotation", "-e",
                os.path.join(self.artifacts_dir, "exLOVD_BRCA1_error_variants.txt"),
                "-s", "exLOVD"]

        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(ConvertEXLOVDBRCA1ExtractToVCF)
class ConvertEXLOVDBRCA2ExtractToVCF(DefaultPipelineTask):

    def output(self):
        return luigi.LocalTarget(self.ex_lovd_file_dir + "/exLOVD_brca2.hg19.vcf")

    def run(self):
        os.chdir(lovd_method_dir)

        args = ["./lovd2vcf.py", "-i", self.ex_lovd_file_dir + "/BRCA2.txt", "-o",
                self.output().path, "-a",
                "exLOVDAnnotation", "-e",
                os.path.join(self.artifacts_dir, "exLOVD_BRCA2_error_variants.txt"),
                "-s", "exLOVD"]

        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(ConvertEXLOVDBRCA2ExtractToVCF)
class ConcatenateEXLOVDVCFFiles(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(self.ex_lovd_file_dir + "/exLOVD_brca12.hg19.vcf")

    def run(self):
        args = ["vcf-concat", self.ex_lovd_file_dir + "/exLOVD_brca1.hg19.vcf",
                self.ex_lovd_file_dir + "/exLOVD_brca2.hg19.vcf"]

        pipeline_utils.run_process(args, redirect_stdout_path=self.output().path)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(ConcatenateEXLOVDVCFFiles)
class CrossmapConcatenatedEXLOVDData(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(self.ex_lovd_file_dir + "/exLOVD_brca12.hg38.vcf")

    def run(self):
        brca_resources_dir = self.cfg.resources_dir

        args = ["CrossMap.py", "vcf",
                brca_resources_dir + "/hg19ToHg38.over.chain.gz",
                self.input().path,
                brca_resources_dir + "/hg38.fa",
                self.output().path]

        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(CrossmapConcatenatedEXLOVDData)
class SortEXLOVDOutput(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(
            self.ex_lovd_file_dir + "/exLOVD_brca12.sorted.hg38.vcf")

    def run(self):
        args = ["vcf-sort", self.input().path]
        pipeline_utils.run_process(args, redirect_stdout_path=self.output().path)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(SortEXLOVDOutput)
class CopyEXLOVDOutputToOutputDir(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(
            self.cfg.output_dir + "/exLOVD_brca12.sorted.hg38.vcf")

    def run(self):
        copy(self.input().path, self.cfg.output_dir)
        pipeline_utils.check_file_for_contents(self.output().path)


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
                                     default='https://databases.lovd.nl/shared/export/',
                                     description='URL to download shared LOVD BRCA data from')


    def output(self):
        if len(str(self.lovd_data_file)) != 0:
            return luigi.LocalTarget(str(luigi.LocalTarget(str(self.lovd_data_file))))
        else:
            output = {}

            for symbol in pipeline_utils.get_lovd_symbols(self.cfg.gene_metadata['symbol']):
                output[symbol] = luigi.LocalTarget(f"{self.lovd_file_dir}/{symbol}.txt")

            return output

    def run(self):
        for symbol in pipeline_utils.get_lovd_symbols(self.cfg.gene_metadata['symbol']):
            pipeline_utils.create_path_if_nonexistent(
                os.path.dirname(self.output()[symbol].path))
            data = pipeline_utils.urlopen_with_retry(
                self.shared_lovd_data_url + symbol).read()
            with open(self.output()[symbol].path, "wb") as f:
                f.write(data)


@requires(DownloadLOVDInputFile)
class ConcatenateLOVDData(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(f"{self.lovd_file_dir}/lovd_concatenated_genes.txt")

    def run(self):
        pipeline_utils.concatenate_files_with_identical_header_rows(self.lovd_file_dir, self.output().path)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(ConcatenateLOVDData)
class NormalizeLOVDSubmissions(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(f"{self.lovd_file_dir}/LOVD_normalized.tsv")

    def run(self):
        os.chdir(lovd_method_dir)
        args = ["python", "normalizeLOVDSubmissions.py", "-i",
                self.input().path, "-o",
                self.output().path, "-b",
                f"{self.artifacts_dir}/LOVD_bracket_variants.tsv"]
        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(NormalizeLOVDSubmissions)
class CombineEquivalentLOVDSubmissions(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(f"{self.lovd_file_dir}/LOVD_normalized_combined.tsv")

    def run(self):
        os.chdir(lovd_method_dir)
        args = ["python", "combineEquivalentVariantSubmissions.py", "-i",
                        self.input().path, "-o",
                        self.output().path]
        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(CombineEquivalentLOVDSubmissions)
class ConvertSharedLOVDToVCF(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(f"{self.lovd_file_dir}/sharedLOVD.hg19.vcf")

    def run(self):
        os.chdir(lovd_method_dir)
        args = ["python", "lovd2vcf.py", "-i", self.input().path, "-o",
                self.output().path, "-a", "sharedLOVDAnnotation", "-e",
                f"{self.artifacts_dir}/LOVD_error_variants.txt",
                "-s", "LOVD"]
        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(ConvertSharedLOVDToVCF)
class CrossmapConcatenatedSharedLOVDData(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(f"{self.lovd_file_dir}/sharedLOVD.hg38.vcf")

    def run(self):
        brca_resources_dir = self.cfg.resources_dir

        args = ["CrossMap.py", "vcf",
                brca_resources_dir + "/hg19ToHg38.over.chain.gz",
                self.input().path,
                brca_resources_dir + "/hg38.fa",
                self.output().path]
        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(CrossmapConcatenatedSharedLOVDData)
class SortSharedLOVDOutput(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(f"{self.lovd_file_dir}/sharedLOVD.sorted.hg38.vcf")

    def run(self):
        args = ["vcf-sort", self.input().path]
        pipeline_utils.run_process(args, redirect_stdout_path=self.output().path)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(SortSharedLOVDOutput)
class CopySharedLOVDOutputToOutputDir(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(f"{self.cfg.output_dir}/sharedLOVD.sorted.hg38.vcf")

    def run(self):
        copy(self.input().path, self.cfg.output_dir)
        pipeline_utils.check_file_for_contents(self.output().path)


###############################################
#                    G1K                      #
###############################################


class DownloadStaticG1KData(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(f"{self.g1k_file_dir}/1000G.sorted.hg38.vcf")

    def run(self):
        os.chdir(self.g1k_file_dir)

        g1k_vcf_url = "https://brcaexchange.org/backend/downloads/1000G.sorted.hg38.vcf"
        pipeline_utils.download_file_and_display_progress(g1k_vcf_url)


@requires(DownloadStaticG1KData)
class CopyG1KOutputToOutputDir(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(f"{self.cfg.output_dir}/1000G.sorted.hg38.vcf")

    def run(self):
        copy(self.input().path, self.cfg.output_dir)
        pipeline_utils.check_file_for_contents(self.output().path)


###############################################
#                    EXAC                     #
###############################################


class DownloadStaticExACData(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(
            self.exac_file_dir + "/exac.brca12.sorted.hg38.vcf")

    def run(self):
        os.chdir(self.exac_file_dir)

        exac_vcf_gz_url = "https://brcaexchange.org/backend/downloads/exac.brca12.sorted.hg38.vcf"
        pipeline_utils.download_file_and_display_progress(exac_vcf_gz_url)


@requires(DownloadStaticExACData)
class CopyEXACOutputToOutputDir(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(
            self.cfg.output_dir + "/exac.brca12.sorted.hg38.vcf")

    def run(self):
        copy(self.input().path,
             self.cfg.output_dir)

        pipeline_utils.check_file_for_contents(self.output().path)


###############################################
#                  ENIGMA                     #
###############################################

@requires(ConvertLatestClinvarDataToXML)
class FilterEnigmaAssertions(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.enigma_file_dir, 'enigma_clinvar.xml'))

    def run(self):
        os.chdir(clinvar_method_dir)

        args = ["python", "filter_enigma_data.py", self.input().path,
                self.output().path]

        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(FilterEnigmaAssertions)
class ExtractEnigmaFromClinvar(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.enigma_file_dir, 'enigma_from_clinvar.tsv'))

    def run(self):
        os.chdir(clinvar_method_dir)

        genes_opts = [ s for g in self.cfg.gene_metadata['symbol'] for s in ['--gene', g]]

        args = ["python", "enigma_from_clinvar.py", self.input().path,
                self.output().path,
                '--logs', os.path.join(self.enigma_file_dir, 'enigma_from_clinvar.log')
               ] + genes_opts

        pipeline_utils.run_process(args)


@requires(ExtractEnigmaFromClinvar)
class CopyEnigmaOutputToOutputDir(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(self.cfg.output_dir + "/enigma_from_clinvar.tsv")

    def run(self):
        copy(self.input().path, self.cfg.output_dir)
        pipeline_utils.check_file_for_contents(self.output().path)


###############################################
#             FUNCTIONAL ASSAYS               #
###############################################


class DownloadFunctionalAssaysInputFile(DefaultPipelineTask):
    findlay_BRCA1_ring_function_scores_url = luigi.Parameter(default='https://brcaexchange.org/backend/downloads/ENIGMA_BRCA12_FunctionalAssays_latest_with_function_scores.tsv',
                                            description='URL to download functional assay data from')

    def output(self):
        return luigi.LocalTarget(self.assays_dir + "/ENIGMA_BRCA12_FunctionalAssays_latest_with_function_scores.tsv")

    def run(self):
        data = pipeline_utils.urlopen_with_retry(self.findlay_BRCA1_ring_function_scores_url).read()
        with open(self.output().path, "wb") as f:
            f.write(data)


@requires(DownloadFunctionalAssaysInputFile)
class ConvertFunctionalAssaysToVCF(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(self.assays_dir + "/ENIGMA_BRCA12_functional_assays_scores.hg19.vcf")

    def run(self):
        os.chdir(functional_assays_method_dir)

        args = ["python", "convert_functional_assay_tsv_to_vcf.py", "-v", "-i", self.input().path, "-o",
                self.output().path, "-a", "functionalAssayAnnotation",
                "-l", self.artifacts_dir + "/functional_assays_error_variants.log",
                "-s", "ENIGMABRCA12FunctionalAssaysFunctionScores"]

        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(ConvertFunctionalAssaysToVCF)
class CrossmapFunctionalAssays(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.assays_dir, "ENIGMA_BRCA12_functional_assays_scores.hg38.vcf"))

    def run(self):
        brca_resources_dir = self.cfg.resources_dir

        args = ["CrossMap.py", "vcf", brca_resources_dir + "/hg19ToHg38.over.chain.gz",
                self.input().path, brca_resources_dir + "/hg38.fa",
                self.output().path]

        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(CrossmapFunctionalAssays)
class SortFunctionalAssays(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.assays_dir, "ENIGMA_BRCA12_functional_assays_scores.sorted.hg38.vcf"))

    def run(self):
        args = ["vcf-sort", self.input().path]

        pipeline_utils.run_process(args, redirect_stdout_path=self.output().path)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(SortFunctionalAssays)
class CopyFunctionalAssaysOutputToOutputDir(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.cfg.output_dir, "ENIGMA_BRCA12_functional_assays_scores.sorted.hg38.vcf"))

    def run(self):
        copy(self.input().path, self.cfg.output_dir)
        pipeline_utils.check_file_for_contents(self.output().path)


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
        yield CopyFunctionalAssaysOutputToOutputDir()

    def output(self):
        return {'merged': luigi.LocalTarget(os.path.join(self.artifacts_dir, "merged.tsv")),
                'reports': luigi.LocalTarget(os.path.join(self.artifacts_dir, "reports.tsv"))}

    def run(self):
        os.chdir(data_merging_method_dir)

        args = ["python", "variant_merging.py", "-i", self.cfg.output_dir + "/",
                "-o",
                self.artifacts_dir + '/', "-a",
                self.artifacts_dir + '/', "-v",
                "-c", self.cfg.gene_config_path]

        pipeline_utils.run_process(args)

        pipeline_utils.check_file_for_contents(self.output()['merged'].path)
        pipeline_utils.check_file_for_contents(self.output()['reports'].path)


@requires(MergeVCFsIntoTSVFile)
class AggregateMergedOutput(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.artifacts_dir, "aggregated.tsv"))

    def run(self):
        os.chdir(data_merging_method_dir)

        args = ["python", "aggregate_across_columns.py",
                "-i", self.input()['merged'].path,
                "-o", self.output().path]

        pipeline_utils.run_process(args)

        pipeline_utils.check_input_and_output_tsvs_for_same_number_variants(
            self.input()['merged'].path,
            self.output().path)


@requires(AggregateMergedOutput)
class BuildAggregatedOutput(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.artifacts_dir, "built.tsv"))

    def run(self):
        brca_resources_dir = self.cfg.resources_dir
        os.chdir(data_merging_method_dir)

        args = ["python", "brca_pseudonym_generator.py",
                self.input().path,
                self.output().path,
                "--log-path", os.path.join(self.artifacts_dir, "brca-pseudonym-generator.log"),
                "--config-file", self.cfg.gene_config_path,
                "--resources", brca_resources_dir]

        pipeline_utils.run_process(args)

        pipeline_utils.check_input_and_output_tsvs_for_same_number_variants(
            self.input().path,
            self.output().path)


@requires(BuildAggregatedOutput)
class AppendCAID(DefaultPipelineTask):

    def output(self):
        artifacts_dir = self.cfg.output_dir + "/release/artifacts/"
        return luigi.LocalTarget(artifacts_dir + "built_with_ca_ids.tsv")

    def run(self):
        release_dir = self.cfg.output_dir + "/release/"
        artifacts_dir = release_dir + "artifacts/"
        brca_resources_dir = self.cfg.resources_dir
        os.chdir(data_merging_method_dir)

        args = ["python", "get_ca_id.py", "-i",
                artifacts_dir + "built.tsv", "-o",
                artifacts_dir + "/built_with_ca_ids.tsv", "-l",
                artifacts_dir + "/get_ca_id.log"]
        print("Running get_ca_id.py with the following args: %s" % (
            args))
        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        pipeline_utils.print_subprocess_output_and_error(sp)

        pipeline_utils.check_input_and_output_tsvs_for_same_number_variants(
            artifacts_dir + "built.tsv",
            artifacts_dir + "built_with_ca_ids.tsv")


@requires(AppendCAID)
class AppendMupitStructure(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.artifacts_dir, "built_with_mupit.tsv"))

    def run(self):
        os.chdir(data_merging_method_dir)

        args = ["python", "getMupitStructure.py", "-i",
                self.input().path, "-o",
                self.output().path]

        pipeline_utils.run_process(args)

        pipeline_utils.check_input_and_output_tsvs_for_same_number_variants(
            self.input().path,
            self.output().path)


@requires(AppendMupitStructure)
class CalculatePriors(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.artifacts_dir, "built_with_priors.tsv"))

    def run(self):
        artifacts_dir_host = self.cfg.output_dir_host + "/release/artifacts/"
        os.chdir(priors_method_dir)

        args = ['bash', 'calcpriors.sh', self.cfg.priors_references_dir,
                artifacts_dir_host, 'built_with_mupit.tsv',
                'built_with_priors.tsv', self.cfg.priors_docker_image_name]

        pipeline_utils.run_process(args)

        pipeline_utils.check_input_and_output_tsvs_for_same_number_variants(
            self.input().path,
            self.output().path)


@requires(CalculatePriors)
class FilterBlacklistedPriors(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.artifacts_dir, "built_with_priors_clean.tsv"))

    def run(self):
        os.chdir(priors_filter_method_dir)

        args = ["python", "filterBlacklistedVars.py",
                "--output", self.output().path,
                "--blacklisted_vars", "blacklisted_vars.txt",
                "filter",
                self.input().path]

        pipeline_utils.run_process(args)

        # we only clear a few columns; we shouldn't be gaining or losing any variants
        pipeline_utils.check_input_and_output_tsvs_for_same_number_variants(
            self.input().path,
            self.output().path)


@requires(FilterBlacklistedPriors)
class AppendVRId(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.artifacts_dir, "built_with_vr_ids.tsv"))

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

        pipeline_utils.run_process(args)

        # we shouldn't be gaining or losing any variants
        pipeline_utils.check_input_and_output_tsvs_for_same_number_variants(
            self.input().path,
            self.output().path)


@requires(bayesdel_processing.AddBayesdelScores)
class FindMissingReports(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.artifacts_dir, "missing_reports.log"))

    def run(self):
        os.chdir(data_merging_method_dir)

        args = ["python", "check_for_missing_reports.py", "-b",
                self.input().path, "-r",
                self.artifacts_dir,
                "-a", self.artifacts_dir, "-v"]

        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(bayesdel_processing.AddBayesdelScores)
class RunDiffAndAppendChangeTypesToOutput(DefaultPipelineTask):
    def _extract_release_date(self, version_json):
        with open(version_json, 'r') as f:
            j = json.load(f)
            return datetime.datetime.strptime(j['date'], '%Y-%m-%d')

    def output(self):
        return {'built_with_change_types': luigi.LocalTarget(
            os.path.join(self.release_dir, "built_with_change_types.tsv")),
                'removed': luigi.LocalTarget(os.path.join(self.diff_dir , "removed.tsv")),
                'added': luigi.LocalTarget(os.path.join(self.diff_dir , "added.tsv")),
                'added_data': luigi.LocalTarget(os.path.join(self.diff_dir , "added_data.tsv")),
                'diff': luigi.LocalTarget(os.path.join(self.diff_dir , "diff.txt")),
                'diff_json': luigi.LocalTarget(os.path.join(self.diff_dir , "diff.json")),
                'README': luigi.LocalTarget(os.path.join(self.diff_dir , "README.txt"))}

    def run(self):
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
                self.input().path, "--v1",
                previous_data_path,
                "--removed", os.path.join(self.diff_dir, "removed.tsv"), "--added",
                os.path.join(self.diff_dir, "added.tsv"), "--added_data",
                os.path.join(self.diff_dir, "added_data.tsv"), "--diff", os.path.join(self.diff_dir, "diff.txt"),
                "--diff_json", os.path.join(self.diff_dir, "diff.json"),
                "--output", os.path.join(self.release_dir, "built_with_change_types.tsv"),
                "--artifacts_dir", self.artifacts_dir,
                "--diff_dir", self.diff_dir, "--v1_release_date",
                previous_release_date_str, "--reports", "False"]

        pipeline_utils.run_process(args)

        shutil.rmtree(tmp_dir)  # cleaning up

        pipeline_utils.check_input_and_output_tsvs_for_same_number_variants(
            os.path.join(self.artifacts_dir, "built_with_vr_ids.tsv"),
            os.path.join(self.release_dir, "built_with_change_types.tsv"))


@requires(RunDiffAndAppendChangeTypesToOutput)
class RunDiffAndAppendChangeTypesToOutputReports(DefaultPipelineTask):
    def _extract_release_date(self, version_json):
        with open(version_json, 'r') as f:
            j = json.load(f)
            return datetime.datetime.strptime(j['date'], '%Y-%m-%d')

    def output(self):
        return {'reports_with_change_types': luigi.LocalTarget(
            os.path.join(self.release_dir, "reports_with_change_types.tsv")),
                'removed_reports': luigi.LocalTarget(
                    os.path.join(self.diff_dir, "removed_reports.tsv")),
                'added_reports': luigi.LocalTarget(
                    os.path.join(self.diff_dir, "added_reports.tsv")),
                'added_data_reports': luigi.LocalTarget(
                    os.path.join(self.diff_dir, "added_data_reports.tsv")),
                'diff_reports': luigi.LocalTarget(
                    os.path.join(self.diff_dir, "diff_reports.txt")),
                'diff_json_reports': luigi.LocalTarget(
                    os.path.join(self.diff_dir, "diff_reports.json")),
                'README': luigi.LocalTarget(os.path.join(self.diff_dir, "README.txt"))}

    def run(self):
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
                os.path.join(self.artifacts_dir, "reports.tsv"), "--v1", previous_data_path,
                "--removed", self.output()['removed_reports'].path, "--added",
                self.output()['added_reports'].path, "--added_data",
                self.output()['added_data_reports'].path, "--diff",
                self.output()['diff_reports'].path, "--diff_json",
                self.output()['diff_json_reports'].path,
                "--output", self.output()['reports_with_change_types'].path,
                "--artifacts_dir", self.artifacts_dir,
                "--diff_dir", self.diff_dir, "--v1_release_date",
                previous_release_date_str, "--reports", "True"]

        pipeline_utils.run_process(args)

        shutil.rmtree(tmp_dir)  # cleaning up

        pipeline_utils.check_input_and_output_tsvs_for_same_number_variants(
            os.path.join(self.artifacts_dir, "reports.tsv"),
            os.path.join(self.release_dir, "reports_with_change_types.tsv"))


@requires(RunDiffAndAppendChangeTypesToOutput)
class LinkBuiltFinal(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.release_dir, "built_final.tsv"))

    def run(self):
        input_path = self.input()['built_with_change_types'].path
        # create relative symlink to have a permanent pointer to what currently the final output of the pipeline is
        relative_input = os.path.relpath(input_path, os.path.dirname(self.output().path))
        os.symlink(relative_input, self.output().path)


@requires(LinkBuiltFinal)
class GenerateVariantsOutputFile(DefaultPipelineTask):
    VAR_OUTPUT_FILE_KEY = 'var_output_file'
    VAR_OUTPUT_METADATA_FILE_KEY = 'var_output_metadata_file'

    def output(self):
        return {
            self.VAR_OUTPUT_FILE_KEY: luigi.LocalTarget(os.path.join(self.cfg.output_dir, "variants_output.tsv")),
            self.VAR_OUTPUT_METADATA_FILE_KEY: luigi.LocalTarget(os.path.join(self.cfg.output_dir, "variants_output_field_metadata.tsv")),
        }

    def run(self):
        os.chdir(data_merging_method_dir)

        in_file = self.input().path
        args = ["python", "generate_variants_output_file.py",
                in_file,
                field_metadata_path,
                field_metadata_path_additional,
                self.output()[self.VAR_OUTPUT_FILE_KEY].path,
                self.output()[self.VAR_OUTPUT_METADATA_FILE_KEY].path
                ]

        pipeline_utils.run_process(args)

        pipeline_utils.check_input_and_output_tsvs_for_same_number_variants(
            in_file,
            self.output()[self.VAR_OUTPUT_FILE_KEY].path)


@requires(GenerateVariantsOutputFile)
class GenerateReleaseNotes(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.metadata_dir, "version.json"))

    def run(self):
        os.chdir(data_merging_method_dir)

        args = ["python", "buildVersionMetadata.py", "--date",
                str(self.cfg.date), "--notes", self.cfg.release_notes,
                "--output", self.output().path]

        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(GenerateReleaseNotes)
class TopLevelReadme(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.cfg.output_dir, "README.txt"))

    def run(self):
        top_level_readme_src = os.path.abspath(
            os.path.join(os.path.realpath(__file__), os.pardir, os.pardir,
                         "top_level_readme.txt"))

        shutil.copyfile(top_level_readme_src, self.output().path)


@requires(TopLevelReadme)
class DataDictionary(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.release_dir, "field_metadata.tsv"))

    def run(self):
        data_dictionary_src = os.path.abspath(
            os.path.join(os.path.realpath(__file__), os.pardir, os.pardir, "field_metadata.tsv"))

        metadata = utils.read_tsv_as_dataframe(data_dictionary_src)
        metadata_public = (metadata.drop(columns=[generate_variants_output_file.VARIANTS_OUTPUT_MAPPING_FIELD])
                                   .sort_values(generate_variants_output_file.FIELD_NAME_FIELD))

        utils.write_dataframe_as_tsv(metadata_public, self.output().path)


@requires(DataDictionary, FindMissingReports, RunDiffAndAppendChangeTypesToOutputReports)
class GenerateMD5Sums(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.cfg.output_dir, "md5sums.txt"))

    def run(self):
        os.chdir(utilities_method_dir)

        workflow_dir = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir))

        keep_list_file_path = os.path.join(workflow_dir, "tarball_files_keep_list.txt")
        discard_list_file_path = os.path.join(workflow_dir, "tarball_files_discard_list.txt")
        args = ["python", "generateMD5Sums.py",
                "-i", self.cfg.output_dir,
                "-o", self.output().path,
                "--keepListFilePath", keep_list_file_path,
                "--discardListFilePath", discard_list_file_path]

        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(GenerateMD5Sums)
class GenerateReleaseArchive(DefaultPipelineTask):
    def output(self):
        # Format archive filename as release-mm-dd-yy.tar.gz
        archive_name = f'release-{self.cfg.date.strftime("%x").replace("/", "-")}.tar.gz'
        return luigi.LocalTarget(Path(self.cfg.output_dir).parent / archive_name)

    def run(self):
        # parse md5sum list.
        with open(self.input().path) as f:
            file_list = [line.strip().split(' ')[-1] for line in f.readlines()]

        # include md5sum in tar as well
        file_list.append(Path(self.input().path).name)

        output_dir = Path(self.cfg.output_dir)
        with tarfile.open(self.output().path, "w:gz") as tar:
            for file in file_list:
                file_path = output_dir / file
                tar.add(file_path, arcname=Path(output_dir.name) / file)

