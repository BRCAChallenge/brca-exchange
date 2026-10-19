import os
import sys

import luigi
from luigi.util import requires

luigi.auto_namespace(scope=__name__)

from workflow.pipeline_common import PipelineParams
from workflow import pipeline_utils
from workflow.variant_assembly import VCFAssembly, VCFAssemblyTask, DownloadGnomADCoverage

_pipeline_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


###############################################
#           RUN VEP ON ALL VARIANTS           #
###############################################

@requires(VCFAssembly)
class AnalyzeVEP(VCFAssemblyTask):
    """Run VEP on all variants and populate analysis_vep (variant_class, variant_type,
    consequences, hgvsp, ptc_genomic_pos)."""

    vep_server_url = luigi.Parameter(
        default='http://localhost:8888',
        description='Base URL of the local VEP REST server')

    def output(self):
        return luigi.LocalTarget(os.path.join(self.vcf_dir, 'analyze_vep.done'))

    def run(self):
        script = os.path.join(_pipeline_dir, 'variant_analysis', 'run_vep_analysis.py')
        args = [sys.executable, script, '--vep-url', self.vep_server_url, '--schema', self.cfg.db_schema]
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
        args = [sys.executable, script, '--victor-vcf', victor_vcf.path, '--schema', self.cfg.db_schema]
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
        args = [sys.executable, script, '--output', self.output().path, '--schema', self.cfg.db_schema]
        self._run_process_with_pipeline_path(args)


@requires(ExportVariantsToVCF)
class GenerateSpliceAIScores(VCFAssemblyTask):
    """Run SpliceAI on all unscored variants and merge with previous scores."""

    genome_fa = luigi.Parameter(
        description='Path to hg38.fa reference genome')
    previous_spliceai_vcf = luigi.OptionalParameter(
        default=None,
        description='Path to SpliceAI-scored VCF from the previous release; '
                     'omit for a from-scratch run to score all variants')
    spliceai_batch_size = luigi.IntParameter(
        default=1000,
        description='Max variants per SpliceAI batch')
    spliceai_depth = luigi.IntParameter(
        default=4999,
        description='SpliceAI search depth (-D)')
    spliceai_annotation = luigi.Parameter(
        default=os.path.join(_pipeline_dir, 'insilico', 'spliceai_annotations',
                             'brca_mane_grch38.txt'),
        description='SpliceAI gene annotation (-A); defaults to the MANE '
                    'Select transcripts for BRCA1/2 rather than SpliceAI\'s '
                    'packaged GENCODE V24 table')
    spliceai_precision = luigi.IntParameter(
        default=3,
        description='Decimal places for the SpliceAI delta scores')

    def output(self):
        return luigi.LocalTarget(os.path.join(self.artifacts_dir, 'variants_with_splice_ai.vcf'))

    def run(self):
        import tempfile
        script = os.path.join(_pipeline_dir, 'insilico', 'add_spliceai_scores_for_new_variants.py')
        tmp_dir = tempfile.mkdtemp()
        args = [
            sys.executable, script,
            '-a', self.input().path,
            '-b', str(self.spliceai_batch_size),
            '-d', str(self.spliceai_depth),
            '-f', self.genome_fa,
            '-g', self.spliceai_annotation,
            '-o', self.output().path,
            '-p', str(self.spliceai_precision),
            '-t', tmp_dir,
        ]
        if self.previous_spliceai_vcf:
            args += ['-s', self.previous_spliceai_vcf]
        self._run_process_with_pipeline_path(args)


@requires(GenerateSpliceAIScores)
class AnalyzeSpliceAI(VCFAssemblyTask):
    """Populate analysis_spliceai from the SpliceAI-scored VCF."""

    overwrite = luigi.BoolParameter(
        default=False,
        description='Re-score variants already present in analysis_spliceai')

    def output(self):
        return luigi.LocalTarget(os.path.join(self.vcf_dir, 'analyze_spliceai.done'))

    def run(self):
        script = os.path.join(_pipeline_dir, 'variant_analysis', 'run_spliceai_analysis.py')
        args = [sys.executable, script, '--spliceai-vcf', self.input().path, '--schema', self.cfg.db_schema]
        if self.overwrite:
            args.append('--overwrite')
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
            sys.executable, script,
            '--processes', str(self.priors_processes),
            '--schema', self.cfg.db_schema,
        ]
        self._run_process_with_pipeline_path(args)
        with open(self.output().path, 'w') as f:
            f.write('done\n')


###############################################
#        COMPUTE + LOAD POPFREQ CODES         #
###############################################

_RESOURCES_DIR = os.path.normpath(os.path.join(_pipeline_dir, '..', '..', 'resources'))


class LCRBed(luigi.ExternalTask):
    """Low-complexity region BED file — expected to be present in the resources directory."""

    lcr_path = luigi.Parameter(
        default=os.path.join(_RESOURCES_DIR, 'LCRFromHengHg38.bed'),
        description='Path to the LCR BED file')

    def output(self):
        return luigi.LocalTarget(self.lcr_path)


@requires(DownloadGnomADCoverage)
class CoverageParquet(VCFAssemblyTask):
    """gnomAD v4.1 joint (exome+genome) coverage parquet.

    Requires DownloadGnomADCoverage so that when gnomAD coverage estimation
    is enabled (--DownloadGnomADCoverage-enabled true), this task waits for
    the freshly-generated parquet rather than racing it; when disabled (the
    default), DownloadGnomADCoverage is already complete against the
    existing resources-directory file and this is a no-op pass-through.
    """

    def output(self):
        return self.input()['joint']


@requires(DownloadGnomADCoverage)
class CoverageParquetV4Exome(VCFAssemblyTask):
    """gnomAD v4.1 exome coverage parquet. See CoverageParquet for why this requires
    DownloadGnomADCoverage rather than reading the resources file directly."""

    def output(self):
        return self.input()['exome']


@requires(DownloadGnomADCoverage)
class CoverageParquetV3Genome(VCFAssemblyTask):
    """gnomAD v3.1 genome coverage parquet. See CoverageParquet for why this requires
    DownloadGnomADCoverage rather than reading the resources file directly."""

    def output(self):
        return self.input()['genome']


@requires(VCFAssembly, CoverageParquet, CoverageParquetV4Exome, CoverageParquetV3Genome, LCRBed)
class AnalyzePopfreq(VCFAssemblyTask):
    """Populate analysis_provisional_evidence_codes with population frequency evidence codes."""

    def output(self):
        return luigi.LocalTarget(os.path.join(self.vcf_dir, 'analyze_popfreq.done'))

    def run(self):
        _, cov_v4_joint, cov_v4_exome, cov_v3_genome, lcr_bed = self.input()
        script = os.path.join(_pipeline_dir, 'variant_analysis', 'run_popfreq_analysis.py')
        args = [
            sys.executable, script,
            '--coverage-v4-joint',  cov_v4_joint.path,
            '--coverage-v4-exome',  cov_v4_exome.path,
            '--coverage-v3-genome', cov_v3_genome.path,
            '--bs1-supporting-faf-threshold', '0.00001',
            '--rare-variant-faf-threshold', '0.00001',
            '--small-indel-size-threshold', '50',
            '--allele-count-rare-variant-threshold', '1',
            '--lcr', lcr_bed.path,
            '--missing-faf-suggests-absence',
            '--method-name', 'popfreq_1.3',
            '--overwrite',
            '--schema', self.cfg.db_schema,
        ]
        self._run_process_with_pipeline_path(args)
        with open(self.output().path, 'w') as f:
            f.write('done\n')


@requires(VCFAssembly, CoverageParquet, CoverageParquetV3Genome)
class AnalyzePopfreqLegacy(VCFAssemblyTask):
    """Populate analysis_provisional_evidence_codes using legacy popfreq_1.2 parameters.

    Uses BS1_Supporting/rare-variant FAF threshold of 0.00002, no LCR filtering,
    and all indels absent from gnomAD v4.1 receive No code met (indel).
    Results are written with method_name='popfreq_1.2'.
    """

    def output(self):
        return luigi.LocalTarget(os.path.join(self.vcf_dir, 'analyze_popfreq_legacy.done'))

    def run(self):
        _, coverage_parquet, cov_v3_genome = self.input()
        script = os.path.join(_pipeline_dir, 'variant_analysis', 'run_popfreq_analysis.py')
        args = [
            sys.executable, script,
            '--coverage-v4-joint',  coverage_parquet.path,
            '--coverage-v3-genome', cov_v3_genome.path,
            '--method-name', 'popfreq_1.2',
            '--bs1-supporting-faf-threshold', '0.00002',
            '--rare-variant-faf-threshold', '0.00002',
            '--small-indel-size-threshold', '0',
            '--allele-count-rare-variant-threshold', '0',
            '--no-lcr',
            '--overwrite',
            '--schema', self.cfg.db_schema,
        ]
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
