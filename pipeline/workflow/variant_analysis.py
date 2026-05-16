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
        script = os.path.join(_pipeline_dir, 'variant_analysis', 'run_vep_analysis.py')
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
        script = os.path.join(_pipeline_dir, 'variant_analysis', 'run_bayesdel_analysis.py')
        args = ['python', script, '--victor-vcf', victor_vcf.path]
        self._run_process_with_pipeline_path(args)
        with open(self.output().path, 'w') as f:
            f.write('done\n')


###############################################
#         GENERATE + LOAD SPLICEAI SCORES     #
###############################################

@requires(VCFAssembly)
class ExportVariantsToVCF(VCFAssemblyTask):
    """Export all GRCh38 variants from the DB to a VCF for SpliceAI input."""

    def output(self):
        return luigi.LocalTarget(os.path.join(self.artifacts_dir, 'all_variants.vcf'))

    def run(self):
        script = os.path.join(_pipeline_dir, 'variant_analysis', 'export_variants_to_vcf.py')
        args = ['python', script, '--output', self.output().path]
        self._run_process_with_pipeline_path(args)


@requires(ExportVariantsToVCF)
class GenerateSpliceAIScores(VCFAssemblyTask):
    """Run SpliceAI on all unscored variants and merge with previous scores."""

    genome_fa = luigi.Parameter(
        description='Path to hg38.fa reference genome')
    previous_spliceai_vcf = luigi.Parameter(
        description='Path to SpliceAI-scored VCF from the previous release')
    spliceai_batch_size = luigi.IntParameter(
        default=1000,
        description='Max variants per SpliceAI batch')
    spliceai_depth = luigi.IntParameter(
        default=4999,
        description='SpliceAI search depth (-D)')

    def output(self):
        return luigi.LocalTarget(os.path.join(self.artifacts_dir, 'variants_with_splice_ai.vcf'))

    def run(self):
        import tempfile
        script = os.path.join(_pipeline_dir, 'insilico', 'add_spliceai_scores_for_new_variants.py')
        tmp_dir = tempfile.mkdtemp()
        args = [
            'python', script,
            '-a', self.input().path,
            '-b', str(self.spliceai_batch_size),
            '-d', str(self.spliceai_depth),
            '-f', self.genome_fa,
            '-g', 'grch38',
            '-o', self.output().path,
            '-s', self.previous_spliceai_vcf,
            '-t', tmp_dir,
        ]
        self._run_process_with_pipeline_path(args)


@requires(GenerateSpliceAIScores)
class AnalyzeSpliceAI(VCFAssemblyTask):
    """Populate analysis_spliceai from the SpliceAI-scored VCF."""

    def output(self):
        return luigi.LocalTarget(os.path.join(self.vcf_dir, 'analyze_spliceai.done'))

    def run(self):
        script = os.path.join(_pipeline_dir, 'variant_analysis', 'run_spliceai_analysis.py')
        args = ['python', script, '--spliceai-vcf', self.input().path]
        self._run_process_with_pipeline_path(args)
        with open(self.output().path, 'w') as f:
            f.write('done\n')


###############################################
#         COMPUTE + LOAD PRIORS SCORES        #
###############################################

@requires(AnalyzeVEP)
class AnalyzePriors(VCFAssemblyTask):
    """Populate analysis_priors with splicing prior probabilities from calcVarPriors."""

    priors_processes = luigi.IntParameter(
        default=8,
        description='Number of parallel calcVarPriors workers')

    def output(self):
        return luigi.LocalTarget(os.path.join(self.vcf_dir, 'analyze_priors.done'))

    def run(self):
        script = os.path.join(_pipeline_dir, 'variant_analysis', 'run_priors_analysis.py')
        args = [
            'python', script,
            '--processes', str(self.priors_processes),
        ]
        self._run_process_with_pipeline_path(args)
        with open(self.output().path, 'w') as f:
            f.write('done\n')


###############################################
#        COMPUTE + LOAD POPFREQ CODES         #
###############################################

@requires(VCFAssembly)
class AnalyzePopfreq(VCFAssemblyTask):
    """Populate analysis_provisional_evidence_codes with population frequency evidence codes."""

    data_dir = luigi.Parameter(
        description='Directory containing df_cov_v4.csv (gnomAD v4 coverage)')
    lcr_bed = luigi.Parameter(
        default='',
        description='Optional BED file of low-complexity regions')

    def output(self):
        return luigi.LocalTarget(os.path.join(self.vcf_dir, 'analyze_popfreq.done'))

    def run(self):
        script = os.path.join(_pipeline_dir, 'variant_analysis', 'run_popfreq_analysis.py')
        args = ['python', script, '--data-dir', self.data_dir]
        if self.lcr_bed:
            args += ['--lcr', self.lcr_bed]
        self._run_process_with_pipeline_path(args)
        with open(self.output().path, 'w') as f:
            f.write('done\n')


###############################################
#               TOP-LEVEL TASK                #
###############################################

@requires(AnalyzeVEP, AnalyzeBayesDel, AnalyzeSpliceAI, AnalyzePriors, AnalyzePopfreq)
class VariantAnalysis(VCFAssemblyTask):
    """Top-level variant analysis task."""

    def output(self):
        return luigi.LocalTarget(os.path.join(self.vcf_dir, 'variant_analysis.done'))

    def run(self):
        with open(self.output().path, 'w') as f:
            f.write('done\n')
