import os
from pathlib import Path
from pipeline_common import DefaultPipelineTask
import pipeline_utils
from CompileVCFFiles import AppendVRId, data_merging_method_dir

import subprocess
import luigi
from luigi.util import requires


@requires(AppendVRId)
class ConvertBuiltToVCF(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(Path(self.cfg.output_dir)/'release'/'artifacts'/'bayesdel.vcf')

    def run(self):
        os.chdir(data_merging_method_dir)

        args = ["python", "bayesdel/convert_merged_variants_to_vcf.py", self.input().path, self.output().path]

        print("Running with the following args: %s" % (
            args))
        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)

        pipeline_utils.print_subprocess_output_and_error(sp)


@requires(ConvertBuiltToVCF)
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
        return luigi.LocalTarget(Path(self.cfg.output_dir) / 'release' / 'artifacts' / 'built_with_bayesdel.tsv')

    def run(self):
        os.chdir(data_merging_method_dir)

        args = ["python", "bayesdel/add_bayesdel_scores_to_built_file.py", '--output', self.output().path,
                '--built-tsv', ConvertBuiltToVCF().input().path] + \
                [ p.path for p in self.input()['paths']]

        print("Running with the following args: %s" % (
            args))
        sp = subprocess.Popen(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)

        pipeline_utils.print_subprocess_output_and_error(sp)
