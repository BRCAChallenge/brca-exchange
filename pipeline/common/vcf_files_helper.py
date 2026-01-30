import logging
import re

import hgvs
import pandas as pd

from .hgvs_utils import HgvsWrapper
from .seq_utils import SeqRepoWrapper
from .variant_utils import VCFVariant

# Pattern to match parentheses in tag names (invalid for VCF INFO fields)
INVALID_TAG_PARENTHESES_PATTERN = re.compile(r'[()]')


def cdna_str_to_genomic_var(cdna_hgvs_str, assembly=HgvsWrapper.GRCh38_Assem, hgvs_wrapper = HgvsWrapper.get_instance(),
                            seq_fetcher=SeqRepoWrapper.get_instance()):
    """
    Convert a CDNA HGVS string into a genomic coordinate.
    Normalizes in genomic space.
    Note, that the seq_fetcher needs to be in the same assembly as given by the `assembly` argument.

    :param cdna_hgvs_str: str
    :param assembly: str, assembly version
    :param hgvs_wrapper: HgvsWrapper instance
    :param seq_fetcher: SeqRepoWrapper instance
    :return: VCF variant object
    """
    # assembly_name of seq_fetcher may contain the patch version, like GRCh38.p11
    assert seq_fetcher.assembly_name.startswith(assembly), \
        "seq_fetcher assembly does not correspond to assembly. " \
        "Got " +  seq_fetcher.assembly_name + " and " + assembly

    hgvs_cdna = hgvs_wrapper.hgvs_parser.parse(cdna_hgvs_str)
    if hgvs_cdna:
        hgvs_g = hgvs_wrapper.nm_to_genomic(hgvs_cdna, assembly)
        if hgvs_g:
            hgvs_g_norm = hgvs.normalizer.Normalizer(hgvs_wrapper.hgvs_dp, shuffle_direction=5).normalize(hgvs_g)
            if hgvs_g_norm:
                return VCFVariant.from_hgvs_obj(hgvs_g_norm, seq_fetcher)
    return None


def normalize_field_value(field_value):
    if not _is_empty_field(field_value):
        if field_value[0] == ';':
            field_value = field_value[1:]
        if field_value[-1] == ';':
            field_value = field_value[:-1]
        if ';' in field_value:
            field_value = field_value.replace(';', '')
    return field_value


def read_vcf_as_dataframe(path):
    return pd.read_csv(path, sep='\t', comment='#', header=None)


def _is_empty_field(field_value):
    return field_value == '' or field_value is None


def has_invalid_parentheses(tag_name):
    """
    Check if a tag name contains parentheses, which are invalid for VCF INFO fields.

    Args:
        tag_name: The tag name to check

    Returns:
        True if the tag name contains parentheses, False otherwise
    """
    return bool(INVALID_TAG_PARENTHESES_PATTERN.search(tag_name))


def sanitize_hgvs_tag_name(tag_name):
    """
    Sanitize an HGVS tag name by removing parentheses.

    VCF INFO field tag names must only contain alphanumeric characters and underscores.
    HGVS nomenclature often uses parentheses (e.g., p.(Arg345Pro)), which are invalid
    as VCF tag names.

    Args:
        tag_name: The tag name to sanitize

    Returns:
        The sanitized tag name with parentheses removed
    """
    return INVALID_TAG_PARENTHESES_PATTERN.sub('', tag_name)


def validate_vcf_info_tags(vcf_reader, file_path):
    """
    Validate VCF INFO field tag names for invalid characters (parentheses).

    Logs warnings for any tags containing parentheses and returns a mapping
    of invalid tag names to their sanitized versions.

    Args:
        vcf_reader: pysam VariantFile object
        file_path: Path to the VCF file (for logging)

    Returns:
        dict: Mapping of invalid tag names to sanitized tag names
    """
    invalid_tags = {}
    header = vcf_reader.header

    for record in header.info:
        tag_name = record
        if has_invalid_parentheses(tag_name):
            sanitized_name = sanitize_hgvs_tag_name(tag_name)
            invalid_tags[tag_name] = sanitized_name
            logging.warning(
                "Invalid HGVS tag name '%s' in %s contains parentheses. "
                "Sanitized version: '%s'",
                tag_name, file_path, sanitized_name
            )

    return invalid_tags
