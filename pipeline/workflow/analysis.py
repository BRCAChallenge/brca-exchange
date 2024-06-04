import os
import subprocess
from pathlib import Path

import luigi
from luigi.util import requires

import workflow.pipeline_utils as pipeline_utils
from workflow.pipeline_common import DefaultPipelineTask, data_merging_method_dir, analysis_method_dir
from workflow import bayesdel_processing


@requires(bayesdel_processing.AddBayesdelScores)
class DownloadResources(DefaultPipelineTask):
    def output(self):
        return(luigi.LocalTarget(Path(self.cfg.output_dir) / 'release' / 'artifacts' / 'analysis'))

    
    def run(self):
        download_url = "https://brcaexchange.org/backend/downloads/variant_scoring/"
        self.inputfiles = [ "df_cov_v2.parquet", "df_cov_v3.parquet",
                            "brca.gnomAD.2.1.1.hg19.flags.tsv",
                            "brca.gnomAD.3.1.1.hg38.flags.tsv"
                            ]
        self.wdir = self.output().path
        if not os.path.exists(self.wdir):
            os.mkdir(self.wdir)
        os.chdir(self.wdir)
        for this_file in self.inputfiles:
            pathname = download_url + "/" + this_file
            pipeline_utils.download_file_and_display_progress(pathname)
        

@requires(DownloadResources)
class runPopfreqAssessment(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(Path(self.cfg.output_dir)/'release'/'artifacts'/'built_with_popfreq.tsv')

    def run(self):
        os.chdir(analysis_method_dir)
        args = ["python", "popfreq_assessment.py",
                "--input", bayesdel_processing.AddBayesdelScores().output().path,
                "--data_dir", DownloadResources().output().path,
                "--output", self.output().path]
        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)
                 
@requires(runPopfreqAssessment)
class runBioinfoPred(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(Path(self.cfg.output_dir)/'release'/'artifacts'/'built_with_bioinfo_pred.tsv')

    def run(self):
        os.chdir(analysis_method_dir)
        args = ["python", "add_bioinfo_pred.py",
                "--input", runPopfreqAssessment().output().path,
                "--output", self.output().path]
        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)
                 

        

