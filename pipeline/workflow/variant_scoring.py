import os
import subprocess
from pathlib import Path

import luigi
from luigi.util import requires

import workflow.pipeline_utils as pipeline_utils
from workflow.pipeline_common import DefaultPipelineTask, data_merging_method_dir, analysis_method_dir



@requires(bayesdel_processing.AddBayesdelScores)
class DownloadResources(DefaultPipelineTask):
    def output(self):
        wdir = Path(self.cfg.output_dir) / 'release' / 'artifacts' / 'variant_scoring_wdir'
        return{'wdir': luigi.LocalTarget(wdir),
               'covV2': "df_cov_v2.parquet",
               'covV3': "df_cov_v2.parquet",
               'flagsV2': "brca.gnomAD.2.1.1.hg19.flags.tsv",
               'flagsV3': "brca.gnomAD.3.1.1.hg38.flags.tsv"
               }
    
    def run(self):
        wdir = self.output()['wdir'].path
        if not os.path.exists(wdir):
            os.mkdir(wdir)
        os.chdir(wdir)
        variant_scoring_data_url = "https://brcaexchange.org/backend/downloads/variant_scoring/*"
        pipeline_utils.download_file_and_display_progress(brca1_data_url)
        

@requires(DownloadResources)
class runPopfreqAssessment(DefaultPipelineTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.artifacts_dir, 'built_with_popfreq.tsv'))

    def run(self):
        os.chdir(analysis_method_dir)
        args = ["python", "popfreq_assessment.py",
                "--input", "bayesdel_processing.AddBayesdelScores.output()",
                "--data_dir", DownloadResources.output()["wdir"],
                "--output", self.output().path]
        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)
                 

        

