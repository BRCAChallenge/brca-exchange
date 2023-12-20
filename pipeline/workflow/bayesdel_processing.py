import os
import subprocess
import tempfile

from pathlib import Path

import luigi
from luigi.util import requires

import workflow.pipeline_utils as pipeline_utils
from workflow.pipeline_common import DefaultPipelineTask, data_merging_method_dir, splice_ai_method_dir


class ConvertBuiltToVCF(DefaultPipelineTask):
    def requires(self):
        from workflow.CompileVCFFiles import AppendVRId
        return AppendVRId()

    def output(self):
        return luigi.LocalTarget(Path(self.cfg.output_dir)/'release'/'artifacts'/'bayesdel.vcf')

    def run(self):
        os.chdir(data_merging_method_dir)

        args = ["python", "bayesdel/convert_merged_variants_to_vcf.py", self.input().path, self.output().path]

        pipeline_utils.run_process(args)


@requires(ConvertBuiltToVCF)
class GenerateSpliceAIData(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.artifacts_dir, "variants_with_splice_ai.vcf"))

    def run(self):
        brca_resources_dir = self.cfg.resources_dir
        os.chdir(splice_ai_method_dir)
              
        #
        # Extract the spliceAI output VCF from the previous release to identify
        # any new variants.  Score the new variants in batches.  Merge the
        # old and new score sets to generate a score set for the
        # current variants
        tmp_dir = tempfile.mkdtemp()
        previous_vcf_path = pipeline_utils.extract_file(
            self.cfg.previous_release_tar, tmp_dir,
            'output/release/artifacts/variants_with_splice_ai.vcf')
        args = ["python", "add_splcieai_scores_for_new_variants.py",
                "-a", self.input().path, "-b", "1000", "-d", "4999",
                "-f", brca_resources_dir + "/hg38.fa", "-f", "grch38",
                "-o", self.output().path, "-s", previous_vcf_path,
                "-t", tmp_dir]
        pipeline_utils.run_process(args)
        shutil.rmtree(tmp_dir)  
                

@requires(GenerateSpliceAIData)
class AddSpliceAI(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.artifacts_dir, 'built_with_spliceai.tsv'))

    def run(self):
        os.chdir(splice_ai_method_dir)

        args = ["python", "add_splice_scores_to_built_file.py", "--vcf", self.input().path,
                '--built-tsv', ConvertBuiltToVCF().input().path, '--output', self.output().path]

        pipeline_utils.run_process(args)


@requires(AddSpliceAI)
class VictorAnnotations(DefaultPipelineTask):
    def output(self):
        wdir = Path(self.cfg.output_dir) / 'release' / 'artifacts' / 'victor_wdir'

        return {'wdir': luigi.LocalTarget(wdir),
                'paths': [luigi.LocalTarget(wdir / f"output.{c + 1}.qc.vcf.gz") for c in range(0, 20)]
                }

    def run(self):
        os.chdir(data_merging_method_dir)

        wdir = self.output()['wdir'].path
        if not os.path.exists(wdir):
            os.mkdir(wdir)

        log_file = Path(self.cfg.output_dir) / 'release' / 'artifacts' / 'victor_annotation.log'

        # Even though we are running a docker container (victor) from within a docker container (BE pipeline), we need to
        # do the file mapping with respect to the host file system not the file system from the BE pipeline container
        vcf_host = Path(self.cfg.output_dir_host) / 'release' / 'artifacts' / 'bayesdel.vcf'
        wdir_host = Path(self.cfg.output_dir_host) / 'release' / 'artifacts' / 'victor_wdir'

        args = ["bash", "bayesdel/run_annotation_docker.sh",
                str(vcf_host),
                str(self.cfg.victor_data_dir),
                str(wdir_host),
                self.cfg.victor_docker_image_name]

        print("Running with the following args: %s" % (
            args))

        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT)

        proc_out, _ = sp.communicate()
        with open(log_file, 'w') as f:
            f.write(proc_out.decode('utf-8'))

        print(f"Wrote log to {log_file}")



@requires(VictorAnnotations)
class AddBayesdelScores(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.artifacts_dir, 'built_with_bayesdel.tsv'))

    def run(self):
        os.chdir(data_merging_method_dir)

        args = ["python", "bayesdel/add_bayesdel_scores_to_built_file.py", '--output', self.output().path,
                '--built-tsv', AddSpliceAI().output().path] + \
                [ p.path for p in self.input()['paths']]

        pipeline_utils.run_process(args)


