import datetime
import os
from shutil import copy

import luigi

from common import config
from workflow import pipeline_utils

DEFAULT_BRCA_RESOURCES_DIR = (os.path.abspath('../brca/brca-resources'))
DEFAULT_OUTPUT_DIR = (
    os.path.abspath('../brca/pipeline-data/data/pipeline_input'))
DEFAULT_FILE_PARENT_DIR = (os.path.abspath('../brca/pipeline-data/data'))


clinvar_method_dir = os.path.abspath('../clinvar')
lovd_method_dir = os.path.abspath('../lovd')
g1k_method_dir = os.path.abspath('../1000_Genomes')
enigma_method_dir = os.path.abspath('../enigma')
functional_assays_method_dir = os.path.abspath('../functional_assays')
data_merging_method_dir = os.path.abspath('../data_merging')
priors_method_dir = os.path.abspath('../splicing')
priors_filter_method_dir = os.path.abspath('../splicingfilter')
utilities_method_dir = os.path.abspath('../utilities')
vr_method_dir = os.path.abspath('../vr')
splice_ai_method_dir = os.path.abspath('../splice_ai')
analysis_method_dir = os.path.abspath('../analysis')

field_metadata_path = os.path.abspath(os.path.join('..', 'field_metadata.tsv'))
field_metadata_path_additional = os.path.abspath(os.path.join('..', 'field_metadata_additional_fields_variant_output_file.tsv'))

class PipelineParams(luigi.Config):
    date = luigi.DateParameter(default=datetime.date.today())
    u = luigi.Parameter(default="UNKNOWN_USER")
    p = luigi.Parameter(default="UNKNOWN_PASSWORD", significant=False)

    resources_dir = luigi.Parameter(default=DEFAULT_BRCA_RESOURCES_DIR,
                                    description='directory to store brca-resources data')

    output_dir = luigi.Parameter(default=DEFAULT_OUTPUT_DIR,
                                 description='directory to store output files')

    output_dir_host = luigi.Parameter(default=DEFAULT_OUTPUT_DIR,
                                      description='directory to store output files wrt to host file system (needed for setting up volume mapping for running docker inside docker)')

    file_parent_dir = luigi.Parameter(default=DEFAULT_FILE_PARENT_DIR,
                                      description='directory to store all individual task related files')

    previous_release_tar = luigi.Parameter(default=str(None), description='path to previous release tar for diffing versions \
                                       and producing change types for variants')

    priors_references_dir = luigi.Parameter(default=str(None),
                                            description='directory to store priors references data')

    priors_docker_image_name = luigi.Parameter(default=str(None),
                                               description='docker image name for priors calculation')

    vr_docker_image_name = luigi.Parameter(default=str(None),
                                               description='docker image name for VR id gathering')

    release_notes = luigi.Parameter(default=str(None),
                                    description='notes for release, must be a .txt file')

    gene_config_path = luigi.Parameter(default=str(None),
                                       description='gene metadata config file')

    seq_repo_dir = luigi.Parameter(default=str(None),
                                   description='directory of seq_repo')

    victor_docker_image_name = luigi.Parameter(default=str(None),
                                               description='docker image name for running the victor annotation pipeline')

    victor_data_dir = luigi.Parameter(default=str(None),
                                               description='data dir with the required data for victor')

    def run(self):
        pass

    @property
    def gene_metadata(self):
        return config.load_config(str(self.gene_config_path))

    __instance = None

    @staticmethod
    def get_instance():
        if not PipelineParams.__instance:
            PipelineParams.__instance = PipelineParams()
        return PipelineParams.__instance


class DefaultPipelineTask(luigi.Task):
    def __init__(self, *args, **kwargs):
        super(DefaultPipelineTask, self).__init__(*args, **kwargs)

        self.cfg = PipelineParams.get_instance()

        self.artifacts_dir = pipeline_utils.create_path_if_nonexistent(self.cfg.output_dir + "/release/artifacts")
        self.release_dir = pipeline_utils.create_path_if_nonexistent(self.cfg.output_dir + "/release/")
        self.diff_dir = pipeline_utils.create_path_if_nonexistent(self.release_dir + "/diff")
        self.metadata_dir = pipeline_utils.create_path_if_nonexistent(self.cfg.output_dir + "/release/metadata/")
        self.esp_file_dir = pipeline_utils.create_path_if_nonexistent(self.cfg.file_parent_dir + "/ESP")
        self.bic_file_dir = pipeline_utils.create_path_if_nonexistent(self.cfg.file_parent_dir + "/BIC")
        self.gnomad_file_dir = pipeline_utils.create_path_if_nonexistent(self.cfg.file_parent_dir + "/gnomAD")
        self.clinvar_file_dir = pipeline_utils.create_path_if_nonexistent(self.cfg.file_parent_dir + "/ClinVar")
        self.ex_lovd_file_dir = pipeline_utils.create_path_if_nonexistent(self.cfg.file_parent_dir + "/exLOVD")
        self.lovd_file_dir = pipeline_utils.create_path_if_nonexistent(self.cfg.file_parent_dir + "/LOVD")
        self.g1k_file_dir = pipeline_utils.create_path_if_nonexistent(self.cfg.file_parent_dir + '/G1K')
        self.exac_file_dir = pipeline_utils.create_path_if_nonexistent(self.cfg.file_parent_dir + '/exac')
        self.enigma_file_dir = pipeline_utils.create_path_if_nonexistent(self.cfg.file_parent_dir + '/enigma')
        self.assays_dir = pipeline_utils.create_path_if_nonexistent(self.cfg.file_parent_dir + '/functional_assays')

    def on_failure(self, exception):
        # renaming files by prefixing filename with "FAILURE_". This way, on rerunning the pipeline the failed task is
        # automatically run again, but instead of just deleting the file, it can still be inspected.
        def _rename_file(path):
            if os.path.exists(path):
                ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                failed_file_name = f"FAILED_{ts}_{os.path.basename(path)}"
                os.rename(path, os.path.join(os.path.dirname(path), failed_file_name))

        if isinstance(self.output(), luigi.LocalTarget):
            _rename_file(self.output().path)
        else:
            for o in self.output().values():
                _rename_file(o.path)

        return super().on_failure(exception)


class CopyOutputToOutputDir(DefaultPipelineTask):
    out_dir = luigi.Parameter()
    req_task = luigi.TaskParameter()

    def __init__(self, out_dir, req_task):
        super(CopyOutputToOutputDir, self).__init__(out_dir=out_dir,
                                                    req_task=req_task)

    def requires(self):
        yield self.req_task

    def output(self):
        if isinstance(self.input()[0], dict):
            outs = {}
            for k, v in self.input()[0].items():
                outs[k] = luigi.LocalTarget(self.out_dir + "/" + os.path.basename(v.path))
            return outs
        else:
            return luigi.LocalTarget(
                self.out_dir + "/" + os.path.basename(self.input()[0].path))

    def _copy_output_to_output_directory(self, required_task_output_path, output_path):
        pipeline_utils.create_path_if_nonexistent(
            os.path.dirname(output_path))

        copy(required_task_output_path, output_path)
        pipeline_utils.check_file_for_contents(output_path)

    def run(self):
        if isinstance(self.output(), dict):
            for k, v in self.output().items():
                self._copy_output_to_output_directory(self.req_task.output()[k].path, v.path)
        else:
            self._copy_output_to_output_directory(self.req_task.output().path, self.output().path)
