import os

import luigi
from luigi.util import requires

luigi.auto_namespace(scope=__name__)

from workflow.pipeline_common import PipelineParams
from workflow import pipeline_utils
from workflow.variant_assembly import VCFAssembly, VCFAssemblyTask

_pipeline_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


###############################################
#           RUN VEP ON ALL VARIANTS           #
###############################################

@requires(VCFAssembly)
class AnalyzeVEP(VCFAssemblyTask):
    """Run VEP on all variants and populate analysis_vep.variant_class."""

    vep_server_url = luigi.Parameter(
        default='http://localhost:8888',
        description='Base URL of the local VEP REST server')

    def output(self):
        return luigi.LocalTarget(os.path.join(self.vcf_dir, 'analyze_vep.done'))

    def run(self):
        script = os.path.join(_pipeline_dir, 'variant_processing', 'run_vep_analysis.py')
        args = ['python', script, '--vep-url', self.vep_server_url]
        self._run_process_with_pipeline_path(args)
        with open(self.output().path, 'w') as f:
            f.write('done\n')


###############################################
#       DOWNLOAD + LOAD BAYESDEL SCORES       #
###############################################

class DownloadVictorAnnotations(VCFAssemblyTask):
    """Download the Victor-annotated VCF containing BayesDel scores."""

    victor_url = luigi.Parameter(
        default='https://brcaexchange.org/backend/downloads/BRCA.ann.all.vcf',
        description='URL of the Victor-annotated VCF')

    def output(self):
        return luigi.LocalTarget(os.path.join(self.artifacts_dir, 'BRCA.ann.all.vcf'))

    def run(self):
        data = pipeline_utils.urlopen_with_retry(self.victor_url).read()
        with open(self.output().path, 'wb') as f:
            f.write(data)


@requires(VCFAssembly, DownloadVictorAnnotations)
class AnalyzeBayesDel(VCFAssemblyTask):
    """Populate analysis_bayesdel with BayesDel_nsfp33a_noAF scores from the Victor VCF."""

    def output(self):
        return luigi.LocalTarget(os.path.join(self.vcf_dir, 'analyze_bayesdel.done'))

    def run(self):
        _, victor_vcf = self.input()
        script = os.path.join(_pipeline_dir, 'variant_processing', 'run_bayesdel_analysis.py')
        args = ['python', script, '--victor-vcf', victor_vcf.path]
        self._run_process_with_pipeline_path(args)
        with open(self.output().path, 'w') as f:
            f.write('done\n')


###############################################
#               TOP-LEVEL TASK                #
###############################################

@requires(AnalyzeVEP, AnalyzeBayesDel)
class VariantAnalysis(VCFAssemblyTask):
    """Top-level variant analysis task."""

    def output(self):
        return luigi.LocalTarget(os.path.join(self.vcf_dir, 'variant_analysis.done'))

    def run(self):
        with open(self.output().path, 'w') as f:
            f.write('done\n')
