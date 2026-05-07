import datetime
import os
import re
import tempfile

import luigi
from luigi.util import requires

luigi.auto_namespace(scope=__name__)

from workflow import pipeline_utils
from workflow.pipeline_common import PipelineParams

_pipeline_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
clinvar_method_dir           = os.path.join(_pipeline_dir, 'clinvar')
functional_assays_method_dir = os.path.join(_pipeline_dir, 'functional_assays')
data_merging_method_dir      = os.path.join(_pipeline_dir, 'data_merging')

# GRCh38 chromosome lengths — keyed by both bare (13) and prefixed (chr13) forms
_GRCh38_CONTIGS = {
    '1': 248956422, '2': 242193529, '3': 198295559, '4': 190214555,
    '5': 181538259, '6': 170805979, '7': 159345973, '8': 145138636,
    '9': 138394717, '10': 133797422, '11': 135086622, '12': 133275309,
    '13': 114364328, '14': 107043718, '15': 101991189, '16': 90338345,
    '17': 83257441, '18': 80373285, '19': 58617616, '20': 64444167,
    '21': 46709983, '22': 50818468, 'X': 156040895, 'Y': 57227415,
    'MT': 16569,
}
_GRCh38_CONTIGS.update({f'chr{k}': v for k, v in _GRCh38_CONTIGS.items()})


class VCFAssemblyTask(luigi.Task):
    """Base task for vcf_assembly. Creates only the directories this pipeline uses."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cfg = PipelineParams.get_instance()

        file_parent_dir = os.path.abspath(self.cfg.file_parent_dir)
        output_dir      = os.path.abspath(self.cfg.output_dir)
        create = pipeline_utils.create_path_if_nonexistent
        self.artifacts_dir    = create(output_dir      + "/release/artifacts")
        self.clinvar_file_dir = create(file_parent_dir + "/ClinVar")
        self.ex_lovd_file_dir = create(file_parent_dir + "/exLOVD")
        self.lovd_file_dir    = create(file_parent_dir + "/LOVD")
        self.enigma_file_dir  = create(file_parent_dir + "/enigma")
        self.assays_dir       = create(file_parent_dir + "/functional_assays")
        self.gnomad_file_dir  = create(file_parent_dir + "/gnomAD")
        self.vcf_dir          = create(file_parent_dir + "/vcf")

    def on_failure(self, exception):
        def _rename_file(path):
            if os.path.exists(path):
                ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                failed_name = f"FAILED_{ts}_{os.path.basename(path)}"
                os.rename(path, os.path.join(os.path.dirname(path), failed_name))

        outputs = self.output()
        if isinstance(outputs, luigi.LocalTarget):
            _rename_file(outputs.path)
        else:
            for o in outputs.values():
                _rename_file(o.path)
        return super().on_failure(exception)

    def _normalize_vcf_header(self, input_vcf, output_vcf):
        """Copy input_vcf to output_vcf, injecting missing ##contig and ##INFO lines."""
        header, data = [], []
        with open(input_vcf) as fh:
            for line in fh:
                if line.startswith('#'):
                    header.append(line)
                else:
                    data.append(line)

        # Inject missing ##contig lines
        declared_contigs = set()
        for line in header:
            m = re.search(r'^##contig=<ID=([^,>]+)', line)
            if m:
                declared_contigs.add(m.group(1))

        referenced_contigs = {row.split('\t')[0] for row in data if row.strip()}
        missing_contigs = [c for c in (referenced_contigs - declared_contigs)
                           if c in _GRCh38_CONTIGS]

        # Inject missing ##INFO lines
        declared_info = set()
        for line in header:
            m = re.search(r'^##INFO=<ID=([^,>]+)', line)
            if m:
                declared_info.add(m.group(1))

        referenced_info = set()
        for row in data:
            if not row.strip():
                continue
            cols = row.split('\t')
            if len(cols) > 7 and cols[7] not in ('.', ''):
                for field in cols[7].split(';'):
                    key = field.split('=')[0]
                    if key:
                        referenced_info.add(key)
        missing_info = sorted(referenced_info - declared_info)

        if not missing_contigs and not missing_info:
            with open(output_vcf, 'w') as fh:
                fh.writelines(header)
                fh.writelines(data)
            return

        chrom_idx = next(i for i, l in enumerate(header) if l.startswith('#CHROM'))
        inject = []
        if missing_contigs:
            inject += [f'##contig=<ID={c},length={_GRCh38_CONTIGS[c]},assembly=GRCh38>\n'
                       for c in sorted(missing_contigs)]
        if missing_info:
            inject += [f'##INFO=<ID={k},Number=.,Type=String,Description="{k}">\n'
                       for k in missing_info]
        header = header[:chrom_idx] + inject + header[chrom_idx:]
        with open(output_vcf, 'w') as fh:
            fh.writelines(header)
            fh.writelines(data)

    def _run_process_with_pipeline_path(self, args, **kwargs):
        """Run a subprocess with _pipeline_dir prepended to PYTHONPATH."""
        prev = os.environ.get('PYTHONPATH', '')
        os.environ['PYTHONPATH'] = _pipeline_dir + ((':' + prev) if prev else '')
        try:
            pipeline_utils.run_process(args, **kwargs)
        finally:
            os.environ['PYTHONPATH'] = prev

    def _vrs_annotate(self, input_vcf, output_vcf_gz, output_pkl):
        with tempfile.NamedTemporaryFile(suffix='.vcf', delete=False) as tmp:
            tmp_path = tmp.name
        try:
            self._normalize_vcf_header(input_vcf, tmp_path)
            args = ["vrs-annotate", "vcf", tmp_path,
                    "--vcf-out", output_vcf_gz,
                    "--pkl-out", output_pkl,
                    "--dataproxy-uri", f"seqrepo+file://{self.cfg.seq_repo_dir}"]
            pipeline_utils.run_process(args)
            pipeline_utils.check_file_for_contents(output_vcf_gz)
        finally:
            os.unlink(tmp_path)


###############################################
#                   CLINVAR                   #
###############################################

class DownloadLatestClinvarData(VCFAssemblyTask):
    def output(self):
        return luigi.LocalTarget(self.clinvar_file_dir + "/ClinVarVCVRelease_00-latest.xml.gz")

    def run(self):
        os.chdir(self.clinvar_file_dir)
        clinvar_data_url = "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarVCVRelease_00-latest.xml.gz"
        pipeline_utils.download_file_and_display_progress(clinvar_data_url)


@requires(DownloadLatestClinvarData)
class ConvertLatestClinvarDataToXML(VCFAssemblyTask):
    def output(self):
        return luigi.LocalTarget(self.clinvar_file_dir + "/ClinVar.xml")

    def run(self):
        os.chdir(clinvar_method_dir)
        genes_opts = [s for g in self.cfg.gene_metadata['symbol'] for s in ['--gene', g]]
        pipeline_utils.run_process(["python", "filter_clinvar.py",
                                    "-i", self.input().path,
                                    "-o", self.output().path] + genes_opts)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(ConvertLatestClinvarDataToXML)
class ConvertClinvarXMLToTXT(VCFAssemblyTask):
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
class ConvertClinvarTXTToVCF(VCFAssemblyTask):
    def output(self):
        return luigi.LocalTarget(self.clinvar_file_dir + "/ClinVar.vcf")

    def run(self):
        os.chdir(data_merging_method_dir)
        args = ["python", "convert_tsv_to_vcf.py",
                "-i", self.input().path,
                "-o", self.output().path,
                "-s", "ClinVar"]
        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(ConvertClinvarTXTToVCF)
class VRSAnnotateClinVar(VCFAssemblyTask):
    def output(self):
        return {
            "vcf": luigi.LocalTarget(f"{self.vcf_dir}/ClinVar.vcf.gz"),
            "pkl": luigi.LocalTarget(f"{self.vcf_dir}/ClinVar.dicts.pkl"),
        }

    def run(self):
        self._vrs_annotate(self.input().path,
                           self.output()["vcf"].path,
                           self.output()["pkl"].path)


###############################################
#                  exLOVD                     #
###############################################

class DownloadStaticExLOVDVCF(VCFAssemblyTask):
    exlovd_static_vcf_url = luigi.Parameter(
        default='https://brcaexchange.org/backend/downloads/vcf/exLOVD.BRCA12.sorted.hg38.vcf',
        description='URL to download static exLOVD vcf from')

    def output(self):
        return luigi.LocalTarget(f"{self.ex_lovd_file_dir}/exLOVD.BRCA12.sorted.hg38.vcf")

    def run(self):
        data = pipeline_utils.urlopen_with_retry(self.exlovd_static_vcf_url).read()
        with open(self.output().path, "wb") as f:
            f.write(data)


@requires(DownloadStaticExLOVDVCF)
class VRSAnnotateEXLOVD(VCFAssemblyTask):
    def output(self):
        return {
            "vcf": luigi.LocalTarget(f"{self.vcf_dir}/exLOVD.BRCA12.sorted.hg38.vcf.gz"),
            "pkl": luigi.LocalTarget(f"{self.vcf_dir}/exLOVD.BRCA12.sorted.hg38.dicts.pkl"),
        }

    def run(self):
        self._vrs_annotate(self.input().path,
                           self.output()["vcf"].path,
                           self.output()["pkl"].path)


##############################################
#               sharedLOVD                   #
##############################################

class DownloadStaticSharedLOVDVCF(VCFAssemblyTask):
    shared_lovd_static_vcf_url = luigi.Parameter(
        default='https://brcaexchange.org/backend/downloads/vcf/LOVD.sorted.hg38.vcf',
        description='URL to download static sharedLOVD vcf from')

    def output(self):
        return luigi.LocalTarget(f"{self.lovd_file_dir}/LOVD.sorted.hg38.vcf")

    def run(self):
        data = pipeline_utils.urlopen_with_retry(self.shared_lovd_static_vcf_url).read()
        with open(self.output().path, "wb") as f:
            f.write(data)


@requires(DownloadStaticSharedLOVDVCF)
class VRSAnnotateSharedLOVD(VCFAssemblyTask):
    def output(self):
        return {
            "vcf": luigi.LocalTarget(f"{self.vcf_dir}/LOVD.sorted.hg38.vcf.gz"),
            "pkl": luigi.LocalTarget(f"{self.vcf_dir}/LOVD.sorted.hg38.dicts.pkl"),
        }

    def run(self):
        self._vrs_annotate(self.input().path,
                           self.output()["vcf"].path,
                           self.output()["pkl"].path)


###############################################
#                  ENIGMA                     #
###############################################

@requires(ConvertLatestClinvarDataToXML)
class FilterEnigmaAssertions(VCFAssemblyTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.enigma_file_dir, 'enigma_clinvar.xml'))

    def run(self):
        os.chdir(clinvar_method_dir)
        args = ["python", "filter_enigma_data.py",
                self.input().path,
                self.output().path]
        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(FilterEnigmaAssertions)
class ExtractEnigmaFromClinvar(VCFAssemblyTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.enigma_file_dir, 'enigma.tsv'))

    def run(self):
        os.chdir(clinvar_method_dir)
        genes_opts = [s for g in self.cfg.gene_metadata['symbol'] for s in ['--gene', g]]
        args = ["python", "enigma_from_clinvar.py",
                self.input().path,
                self.output().path,
                '--logs', os.path.join(self.enigma_file_dir, 'enigma.log'),
               ] + genes_opts
        self._run_process_with_pipeline_path(args)


@requires(ExtractEnigmaFromClinvar)
class ConvertEnigmaToVCF(VCFAssemblyTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.enigma_file_dir, 'enigma.vcf'))

    def run(self):
        os.chdir(data_merging_method_dir)
        args = ["python", "convert_tsv_to_vcf.py",
                "-i", self.input().path,
                "-o", self.output().path,
                "-s", "ENIGMA"]
        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(ConvertEnigmaToVCF)
class VRSAnnotateEnigma(VCFAssemblyTask):
    def output(self):
        return {
            "vcf": luigi.LocalTarget(f"{self.vcf_dir}/enigma.vcf.gz"),
            "pkl": luigi.LocalTarget(f"{self.vcf_dir}/enigma.dicts.pkl"),
        }

    def run(self):
        self._vrs_annotate(self.input().path,
                           self.output()["vcf"].path,
                           self.output()["pkl"].path)


###############################################
#             FUNCTIONAL ASSAYS               #
###############################################

class DownloadFunctionalAssaysInputFile(VCFAssemblyTask):
    findlay_BRCA1_ring_function_scores_url = luigi.Parameter(
        default='https://brcaexchange.org/backend/downloads/ENIGMA_BRCA12_FunctionalAssays_latest_with_function_scores.tsv',
        description='URL to download functional assay data from')

    def output(self):
        return luigi.LocalTarget(
            self.assays_dir + "/ENIGMA_BRCA12_FunctionalAssays_latest_with_function_scores.tsv")

    def run(self):
        data = pipeline_utils.urlopen_with_retry(
            self.findlay_BRCA1_ring_function_scores_url).read()
        with open(self.output().path, "wb") as f:
            f.write(data)


@requires(DownloadFunctionalAssaysInputFile)
class ConvertFunctionalAssaysToVCF(VCFAssemblyTask):
    def output(self):
        return luigi.LocalTarget(
            self.assays_dir + "/ENIGMA_BRCA12_functional_assays_scores.hg19.vcf")

    def run(self):
        os.chdir(functional_assays_method_dir)
        args = ["python", "convert_functional_assay_tsv_to_vcf.py",
                "-v", "-i", self.input().path,
                "-o", self.output().path,
                "-a", "functionalAssayAnnotation",
                "-l", self.artifacts_dir + "/functional_assays_error_variants.log",
                "-s", "ENIGMABRCA12FunctionalAssaysFunctionScores"]
        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(ConvertFunctionalAssaysToVCF)
class CrossmapFunctionalAssays(VCFAssemblyTask):
    def output(self):
        return luigi.LocalTarget(
            os.path.join(self.assays_dir, "ENIGMA_BRCA12_functional_assays_scores.hg38.vcf"))

    def run(self):
        brca_resources_dir = self.cfg.resources_dir
        args = ["CrossMap.py", "vcf",
                brca_resources_dir + "/hg19ToHg38.over.chain.gz",
                self.input().path,
                brca_resources_dir + "/hg38.fa",
                self.output().path]
        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(CrossmapFunctionalAssays)
class SortFunctionalAssays(VCFAssemblyTask):
    def output(self):
        return luigi.LocalTarget(
            os.path.join(self.assays_dir, "ENIGMA_BRCA12_functional_assays_scores.sorted.hg38.vcf"))

    def run(self):
        args = ["vcf-sort", self.input().path]
        pipeline_utils.run_process(args, redirect_stdout_path=self.output().path)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(SortFunctionalAssays)
class VRSAnnotateFunctionalAssays(VCFAssemblyTask):
    def output(self):
        return {
            "vcf": luigi.LocalTarget(f"{self.vcf_dir}/ENIGMA_BRCA12_functional_assays_scores.sorted.hg38.vcf.gz"),
            "pkl": luigi.LocalTarget(f"{self.vcf_dir}/ENIGMA_BRCA12_functional_assays_scores.sorted.hg38.dicts.pkl"),
        }

    def run(self):
        self._vrs_annotate(self.input().path,
                           self.output()["vcf"].path,
                           self.output()["pkl"].path)


###############################################
#                  gnomAD                     #
###############################################

class GnomADTask(VCFAssemblyTask):
    def _download_file(self, url, output_path):
        data = pipeline_utils.urlopen_with_retry(url).read()
        with open(output_path, "wb") as f:
            f.write(data)


class DownloadStaticGnomADVCF(GnomADTask):
    gnomAD_v2_static_vcf_url = luigi.Parameter(
        default='https://brcaexchange.org/backend/downloads/vcf/gnomADv2.sorted.hg38.vcf',
        description='URL to download static gnomAD v2 vcf from')
    gnomAD_v3_static_data_url = luigi.Parameter(
        default='https://brcaexchange.org/backend/downloads/vcf/gnomADv3.sorted.hg38.vcf',
        description='URL to download static gnomAD v3 vcf from')
    gnomad_v4_static_data_url = luigi.Parameter(
        default='https://brcaexchange.org/backend/downloads/vcf/gnomADv4.sorted.hg38.vcf',
        description='URL to download static gnomAD v4 vcf from')

    def output(self):
        return {
            "v2": luigi.LocalTarget(f"{self.gnomad_file_dir}/gnomADv2.sorted.hg38.vcf"),
            "v3": luigi.LocalTarget(f"{self.gnomad_file_dir}/gnomADv3.sorted.hg38.vcf"),
            "v4": luigi.LocalTarget(f"{self.gnomad_file_dir}/gnomADv4.sorted.hg38.vcf"),
        }

    def run(self):
        self._download_file(self.gnomAD_v2_static_vcf_url, self.output()["v2"].path)
        self._download_file(self.gnomAD_v3_static_data_url, self.output()["v3"].path)
        self._download_file(self.gnomad_v4_static_data_url, self.output()["v4"].path)


@requires(DownloadStaticGnomADVCF)
class VRSAnnotateGnomAD(VCFAssemblyTask):
    def output(self):
        return {
            "v2_vcf": luigi.LocalTarget(f"{self.vcf_dir}/gnomADv2.sorted.hg38.vcf.gz"),
            "v2_pkl": luigi.LocalTarget(f"{self.vcf_dir}/gnomADv2.sorted.hg38.dicts.pkl"),
            "v3_vcf": luigi.LocalTarget(f"{self.vcf_dir}/gnomADv3.sorted.hg38.vcf.gz"),
            "v3_pkl": luigi.LocalTarget(f"{self.vcf_dir}/gnomADv3.sorted.hg38.dicts.pkl"),
            "v4_vcf": luigi.LocalTarget(f"{self.vcf_dir}/gnomADv4.sorted.hg38.vcf.gz"),
            "v4_pkl": luigi.LocalTarget(f"{self.vcf_dir}/gnomADv4.sorted.hg38.dicts.pkl"),
        }

    def run(self):
        for version in ("v2", "v3", "v4"):
            self._vrs_annotate(
                self.input()[version].path,
                self.output()[f"{version}_vcf"].path,
                self.output()[f"{version}_pkl"].path,
            )


###############################################
#           LOAD VCFs TO DATABASE             #
###############################################

@requires(
    VRSAnnotateEnigma,
    VRSAnnotateClinVar,
    VRSAnnotateSharedLOVD,
    VRSAnnotateEXLOVD,
    VRSAnnotateGnomAD,
)
class LoadVCFsToDatabase(VCFAssemblyTask):
    """Call the load_vcf Django management command to load all VCFs into the pipeline DB."""

    django_dir = luigi.Parameter(
        default='/data/new_schema/code/website/django',
        description='Directory containing Django manage.py')

    def output(self):
        return luigi.LocalTarget(os.path.join(self.vcf_dir, "load_vcfs_to_db.done"))

    def run(self):
        enigma_in, clinvar_in, lovd_in, exlovd_in, gnomad_in = self.input()
        args = [
            "python", "manage.py", "load_vcf", "--skip-checks",
            "--enigma-vcf",    enigma_in["vcf"].path,
            "--enigma-pkl",    enigma_in["pkl"].path,
            "--clinvar-vcf",   clinvar_in["vcf"].path,
            "--clinvar-pkl",   clinvar_in["pkl"].path,
            "--lovd-vcf",      lovd_in["vcf"].path,
            "--lovd-pkl",      lovd_in["pkl"].path,
            "--exlovd-vcf",    exlovd_in["vcf"].path,
            "--exlovd-pkl",    exlovd_in["pkl"].path,
            "--gnomad-v2-vcf", gnomad_in["v2_vcf"].path,
            "--gnomad-v2-pkl", gnomad_in["v2_pkl"].path,
            "--gnomad-v3-vcf", gnomad_in["v3_vcf"].path,
            "--gnomad-v3-pkl", gnomad_in["v3_pkl"].path,
            "--gnomad-v4-vcf", gnomad_in["v4_vcf"].path,
            "--gnomad-v4-pkl", gnomad_in["v4_pkl"].path,
        ]
        os.chdir(self.django_dir)
        pipeline_utils.run_process(args)
        with open(self.output().path, "w") as f:
            f.write("done\n")


###############################################
#               TOP-LEVEL TASK                #
###############################################

@requires(
    VRSAnnotateClinVar,
    VRSAnnotateEXLOVD,
    VRSAnnotateSharedLOVD,
    VRSAnnotateEnigma,
    VRSAnnotateFunctionalAssays,
    VRSAnnotateGnomAD,
    LoadVCFsToDatabase,
)
class VCFAssembly(VCFAssemblyTask):
    """Runs all source chains, VRS-annotates each VCF, loads to DB, and writes sentinel."""

    def output(self):
        return luigi.LocalTarget(os.path.join(self.vcf_dir, "vcf_assembly.done"))

    def run(self):
        with open(self.output().path, "w") as f:
            f.write("done\n")
