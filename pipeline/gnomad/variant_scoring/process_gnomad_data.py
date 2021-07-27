import functools
import logging
import operator
import shutil
from pathlib import Path

import glow
from pyspark.conf import SparkConf
from pyspark.sql import SparkSession
from pyspark.sql import functions as sf

import constants
import variant_scoring
from common import config as brca_config


def get_spark_session(cores, mb_per_core, spark_tmp_dir):
    spark_tmp_dir.mkdir(exist_ok=True, parents=True)

    driver_mem = cores * mb_per_core + 2000  # + overhead

    spark_cfg = (SparkConf().set("spark.driver.memory", "{}m".format(driver_mem)).
                 set("spark.executor.memory", "{}m".format(mb_per_core)).
                 set("spark.master", "local[{}]".format(cores)).
                 set("spark.sql.execution.arrow.enabled", str(True)).
                 set("spark.jars.packages", 'io.projectglow:glow-spark3_2.12:1.0.1').
                 set("spark.hadoop.io.compression.codecs", "io.projectglow.sql.util.BGZFCodec").
                 set("spark.local.dir", str(spark_tmp_dir))
                 )

    spark = SparkSession.builder.config(conf=spark_cfg).getOrCreate()
    spark = glow.register(spark)
    return spark


def remove_if_exists(path):
    if path.exists():
        shutil.rmtree(path)


def read_raw_coverage_data(input_path, spark):
    return spark.read.csv(str(input_path), header=True, sep='\t', inferSchema=True)


def read_raw_vcf_data(input_paths, spark):
    # using glow
    dfs = [spark.read.format('vcf').load(str(vcf_input)) for vcf_input in input_paths]

    df = dfs[0]
    for d in dfs[1:]:
        df = df.unionAll(d)

    return df


def boundaries_predicate_variants(boundaries):
    # TODO: verify boundary conditions!! (analoguous to brca pipeline)
    chrom_predicate = functools.reduce(operator.or_,
                                       (((sf.col('contigName') == chrom) & (bound[0] <= sf.col('start')) & (
                                               sf.col('end') <= bound[1])) for (chrom, bound) in boundaries.items()))

    return chrom_predicate


# TODO: unify with variants predicate?
def boundaries_predicate_coverage(boundaries):
    return functools.reduce(operator.or_,
                            (((sf.col('chrom') == chrom) &
                              (bound[0] <= sf.col('pos')) &
                              (sf.col('pos') <= bound[1])) for (chrom, bound) in boundaries.items()))


def read_coverage_data_v2(path, spark):
    df_cov = read_raw_coverage_data(str(path), spark).withColumn('chrom', sf.col('chrom').astype('int'))
    return df_cov


def read_coverage_data_v3(path, spark):
    return (read_raw_coverage_data(str(path), spark).
            withColumnRenamed('median_approx', 'median').
            withColumn('chrom', sf.split('locus', ':').getItem(0)).
            withColumn('pos', sf.split('locus', ':').getItem(1).astype('integer')).
            drop('locus').
            where((sf.col('chrom') != "chrX") & (sf.col('chrom') != "chrY")).
            withColumn('chrom', sf.substring('chrom', len('chr') + 1, 2).astype('int'))
            ).select('chrom', 'pos', 'mean',
                     'median', 'over_1', 'over_5', 'over_10', 'over_15',
                     'over_20', 'over_25', 'over_30', 'over_50', 'over_100')


def prepare_variant_data(vcf_paths, boundaries, spark, additional_cols=tuple(constants.additional_cols)):
    spark_df = (read_raw_vcf_data(vcf_paths, spark).
                withColumn('contigName', sf.col('contigName').astype('int'))
                )

    df_var = (spark_df.
              filter(boundaries_predicate_variants(boundaries)).
              select(constants.vcf_mandatory_cols + list(additional_cols))
              ).toPandas()

    df_var = df_var.rename(columns=lambda c: c.replace('INFO_', ''))

    # TODO: need some more explosions
    # TODO: do later in pipeline?
    df_var = (extract_singleton_array_cols(df_var, ['alternateAlleles', 'popmax'] + constants.faf95_col_names).
              rename(columns={'alternateAlleles ': 'alternateAllele'}))

    df_var = variant_scoring.add_name_col(df_var)

    return df_var


def prepare_coverage_data(coverage_path, coverage_reader, boundaries, spark):
    df_cov = (coverage_reader(coverage_path, spark).
              filter(boundaries_predicate_coverage(boundaries))
              ).toPandas()

    return df_cov


def main():
    # input_dir output_dir spark config
    # TODO
    input_dir = Path('/TODO')
    output_dir = Path('/TODO')
    brca_config_dir = Path('/TODO') / 'workflow' / 'gene_config_brca_only.txt'

    df_var_v2_path = output_dir / 'df_var_v2.parquet'
    df_var_v3_path = output_dir / 'df_var_v3.parquet'

    df_cov_v2_path = output_dir / 'df_cov_v2.parquet'
    df_cov_v3_path = output_dir / 'df_cov_v3.parquet'
    gene_config = brca_config.load_config(brca_config_dir)

    # TODO: streamline!
    boundaries37 = {c: (s, e) for _, (c, s, e) in gene_config[['chr', 'start_hg37', 'end_hg37']].iterrows()}
    boundaries38 = {c: (s, e) for _, (c, s, e) in gene_config[['chr', 'start_hg38', 'end_hg38']].iterrows()}

    with get_spark_session(8, 2048, Path('/tmp')) as spark:
        # TODO: keep guards?
        if not df_var_v2_path.exists():
            logging.info("Processing v2 variants")
            df_var_v2 = prepare_variant_data(
                [input_dir / 'gnomad.exomes.r2.1.1.sites.13.vcf.bgz',
                 input_dir / 'gnomad.exomes.r2.1.1.sites.17.vcf.bgz'],
                boundaries37,
                spark)

            df_var_v2.to_parquet(df_var_v2_path)

        if not df_var_v3_path.exists():
            logging.info("Processing v3 variants")
            df_var_v3 = prepare_variant_data(
                [input_dir / 'gnomad.genomes.v3.1.1.sites.chr13.vcf.bgz',
                 input_dir / 'gnomad.genomes.v3.1.1.sites.chr17.vcf.bgz'],
                boundaries38,
                spark)

            df_var_v3.to_parquet(df_var_v3_path)

        logging.info("Processing v2 coverage summaries")

        if not df_cov_v2_path.exists():
            df_cov_v2 = prepare_coverage_data(input_dir / 'gnomad.exomes.coverage.summary.tsv.bgz',
                                              read_coverage_data_v2, boundaries37, spark)
            df_cov_v2.to_parquet(df_cov_v2_path)

        if not df_cov_v3_path.exists():
            df_cov_v3 = prepare_coverage_data(input_dir / 'gnomad.genomes.r3.0.1.coverage.summary.tsv.bgz',
                                              read_coverage_data_v3, boundaries38, spark)
            df_cov_v3.to_parquet(df_cov_v3_path)


if __name__ == "__main__":
    main()
