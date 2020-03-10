from .hgvs_utils import HgvsWrapper
from .variant_utils import VCFVariant
import hgvs
from .seq_utils import SeqRepoWrapper


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

    hgvs_g = hgvs_wrapper.nm_to_genomic(hgvs_cdna, assembly)

    hgvs_g_norm = hgvs.normalizer.Normalizer(hgvs_wrapper.hgvs_dp, shuffle_direction=5).normalize(hgvs_g)
    return VCFVariant.from_hgvs_obj(hgvs_g_norm, seq_fetcher)


def normalize_field_value(field_value):
    if not _is_empty_field(field_value):
        if field_value[0] == ';':
            field_value = field_value[1:]
        if field_value[-1] == ';':
            field_value = field_value[:-1]
        if ';' in field_value:
            field_value = field_value.replace(';', '')
    return field_value


def _is_empty_field(field_value):
    return field_value == '' or field_value is None
