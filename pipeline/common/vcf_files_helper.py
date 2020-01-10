
from .hgvs_utils import HgvsWrapper
from .variant_utils import VCFVariant
import hgvs

# TODO: remove assembly parameter
def cdna_str_to_genomic_var(cdna_hgvs_str, assembly=HgvsWrapper.GRCh38_Assem):
    hgvs_wrapper = HgvsWrapper.get_instance()

    hgvs_cdna = hgvs_wrapper.hgvs_parser.parse(cdna_hgvs_str)

    hgvs_g = hgvs_wrapper.nm_to_genomic(hgvs_cdna, assembly)

    hgvs_g_norm = hgvs.normalizer.Normalizer(hgvs_wrapper.hgvs_dp, shuffle_direction=5).normalize(hgvs_g)
    return VCFVariant.from_hgvs_obj(hgvs_g_norm)

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
