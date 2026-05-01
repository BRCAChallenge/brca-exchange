import os
import re
import tempfile
import json
import io
from operator import __or__
from django.core import serializers
from django.db import connection
from django.db.models import Q
from django.db.models import Value
from django.db.models.functions import Concat
from django.forms.models import model_to_dict
from django.http import JsonResponse, HttpResponse, HttpResponseBadRequest
from django.views.decorators.gzip import gzip_page
from .models import (
    Variant, DataRelease,
    InSilicoPriors, Variant_in_Paper, Paper, VariantRepresentation
)
from django.views.decorators.http import require_http_methods

from datetime import datetime
from operator import itemgetter
import logging
from functools import reduce

DISALLOWED_SEARCH_CHARS = ['\x00']


def releases(request):
    release_id = request.GET.get('release_id')
    if release_id:
        releases = list(DataRelease.objects.filter(id=release_id).values())
    else:
        releases = list(DataRelease.objects.values())
    latest = DataRelease.objects.order_by('-id')[0].id
    response = JsonResponse({"releases": list(releases), "latest": latest})
    response['Access-Control-Allow-Origin'] = '*'
    return response


def variant_counts(request):
    query = Variant.objects.all()
    total_count = query.count()
    brca1_count = query.filter(Gene_Symbol='BRCA1').count()
    brca2_count = query.filter(Gene_Symbol='BRCA2').count()
    query = query.filter(enigma_reports__isnull=False).distinct()
    enigma_count = query.count()
    enigma_pathogenic_count = query.filter(Pathogenicity='Pathogenic').count()
    enigma_benign_count = query.filter(Pathogenicity__contains='Benign').count()
    enigma_likely_benign_count = query.filter(Pathogenicity__contains='Likely benign').count()
    enigma_likely_pathogenic_count = query.filter(Pathogenicity__contains='Likely pathogenic').count()
    query_brca1 = query.filter(Gene_Symbol='BRCA1')
    brca1_enigma_pathogenic_count = query_brca1.filter(Pathogenicity='Pathogenic').count()
    brca1_enigma_benign_count = query_brca1.filter(Pathogenicity__contains='Benign').count()
    brca1_enigma_likely_benign_count = query_brca1.filter(Pathogenicity__contains='Likely benign').count()
    brca1_enigma_likely_pathogenic_count = query_brca1.filter(Pathogenicity__contains='Likely pathogenic').count()
    query_brca2 = query.filter(Gene_Symbol='BRCA2')
    brca2_enigma_pathogenic_count = query_brca2.filter(Pathogenicity='Pathogenic').count()
    brca2_enigma_benign_count = query_brca2.filter(Pathogenicity__contains='Benign').count()
    brca2_enigma_likely_benign_count = query_brca2.filter(Pathogenicity__contains='Likely benign').count()
    brca2_enigma_likely_pathogenic_count = query_brca2.filter(Pathogenicity__contains='Likely pathogenic').count()
    response = JsonResponse({
        "total": total_count,
        "brca1": {
            "total": brca1_count,
            "pathogenic": brca1_enigma_pathogenic_count,
            "benign": brca1_enigma_benign_count,
            "likelyBenign": brca1_enigma_likely_benign_count,
            "likelyPathogenic": brca1_enigma_likely_pathogenic_count },
        "brca2": {
            "total": brca2_count,
            "pathogenic": brca2_enigma_pathogenic_count,
            "benign": brca2_enigma_benign_count,
            "likelyBenign": brca2_enigma_likely_benign_count,
            "likelyPathogenic": brca2_enigma_likely_pathogenic_count },
        "enigma": enigma_count,
        "enigmaPathogenic": enigma_pathogenic_count,
        "enigmaLikelyPathogenic": enigma_likely_pathogenic_count,
        "enigmaBenign": enigma_benign_count,
        "enigmaLikelyBenign": enigma_likely_benign_count })
    response['Access-Control-Allow-Origin'] = '*'
    return response


def variant(request):
    variant_id = request.GET.get('variant_id')

    variant_id_clean = variant_id.lower().strip()
    variant_id_clean = remove_disallowed_chars(variant_id_clean)

    # Accept search by ClinGen Allele Registry id (CA_ID)
    if variant_id_clean.startswith('ca'):
        query = Variant.objects
        variant = query.filter(Q(CA_ID__icontains=variant_id_clean))

        if len(variant) != 1:
            response = JsonResponse({
                "redirect": True,
                "data": variant_id
            })
            return response
        else:
            variant = variant[0]

    else:
        variant = Variant.objects.get(VRS_Digest=variant_id)

    query = Variant.objects.filter(VRS_Digest=variant.VRS_Digest)

    variant_versions = list(map(variant_to_dict, query))
    response = JsonResponse({"data": variant_versions})
    response['Access-Control-Allow-Origin'] = '*'
    return response


def vrid(request):
    vr_id = request.GET.get('vr_id')

    variant = Variant.objects.get(VRS_Digest=vr_id)

    variant_data = variant_to_dict(variant)

    relevant_keys = {
        'id',
        'Gene_Symbol',
        'Reference_Sequence',
        'HGVS_cDNA',
        'BIC_Nomenclature',
        'HGVS_Protein',
        'Protein_Change',
        'CA_ID',
        'VRS',
    }

    response = JsonResponse({"data": {k: variant_data[k] for k in variant_data if k in relevant_keys}})
    response['Access-Control-Allow-Origin'] = '*'
    return response


def variantreps(request):
    vr_reps = list(
        VariantRepresentation.objects.raw("""
        select VR.id, VR."Genomic_Coordinate_hg38", V.id as Variant_id, VR."Description" from data_variantrepresentation VR
        inner join variant V on V."Genomic_Coordinate_hg38" = VR."Genomic_Coordinate_hg38"
        """)
    )

    response = JsonResponse({
        "count": len(vr_reps),
        "data": list(
            {'id': x.variant_id, 'Genomic_Coordinate_hg38': x.Genomic_Coordinate_hg38, 'vr_rep': x.Description}
            for x in vr_reps
        )
    })
    response['Access-Control-Allow-Origin'] = '*'

    return response


def sitemap(request):
    variants = Variant.objects.values_list('VRS_Digest')
    root_links = [
        'https://brcaexchange.org/',
        'https://brcaexchange.org/factsheet',
        'https://brcaexchange.org/help',
        'https://brcaexchange.org/community',
        'https://brcaexchange.org/variants',
        'https://brcaexchange.org/about/thisSite',
        'https://brcaexchange.org/releases',
    ]
    response = HttpResponse(
        ("\n".join(root_links) + "\n") +
        "\n".join("https://brcaexchange.org/variant/%s" % x[0] for x in variants)
    )
    response['Access-Control-Allow-Origin'] = '*'
    response['Content-Type'] = 'text/plain'
    return response


def variant_papers(request):
    variant_id = request.GET.get('variant_id')
    variant = Variant.objects.get(VRS_Digest=variant_id)
    variantpapers = Variant_in_Paper.objects.select_related('Paper').filter(VRS_Digest=variant).all()
    for variantpaper in variantpapers:
        # year of 0000 means year could not be found during a crawl
        if variantpaper.Paper.Year == '0000':
            variantpaper.Paper.Year = "Unknown"
    variantpapers = [dict(model_to_dict(vp.Paper), **{"mentions": vp.mentions, "points": vp.points}) for vp in variantpapers]
    response = JsonResponse({"data": variantpapers}, safe=False)
    response['Access-Control-Allow-Origin'] = '*'
    return response


def variant_to_dict(variant_object):
    variant_dict = model_to_dict(variant_object)

    return variant_dict


@gzip_page
def index(request):
    order_by = request.GET.get('order_by')
    direction = request.GET.get('direction')
    page_size = int(request.GET.get('page_size', '0'))
    page_num = int(request.GET.get('page_num', '0'))
    search_term = request.GET.get('search_term')
    format = request.GET.get('format')
    include = request.GET.getlist('include')
    exclude = request.GET.getlist('exclude')
    filters = request.GET.getlist('filter')
    filter_values = request.GET.getlist('filterValue')
    column = request.GET.getlist('column')
    deleted_count = 0
    synonyms_count = 0

    query = Variant.objects

    if format == 'csv' or format == 'tsv':
        quotes = '\''
    else:
        quotes = ''

    query = apply_sources(query, include, exclude)

    if filters:
        query = apply_filters(query, filter_values, filters, quotes=quotes)

    if search_term:
        query, synonyms_count = apply_search(query, search_term, quotes=quotes)

    if order_by:
        query = apply_order(query, order_by, direction)

    if format == 'csv' or format == 'tsv':
        cursor = connection.cursor()

        f = io.StringIO()
        query = "COPY ({}) TO STDOUT WITH DELIMITER '{}' CSV HEADER".format(query.query, '\t' if format == 'tsv' else ',')
        # HACK to add quotes around search terms
        query = re.sub(r'LIKE UPPER\((.+?)\)', r"LIKE UPPER('\1')", query)
        cursor.copy_expert(query, f)
        f.seek(0)

        response = HttpResponse(f.read(), content_type='text/csv')
        response['Content-Disposition'] = 'attachment;filename="variants.%s"' % format
        return response

    elif format == 'json':
        count = query.count()
        query = select_page(query, page_size, page_num)
        response = JsonResponse({'count': count, 'deletedCount': deleted_count, 'synonyms': synonyms_count, 'data': list(query.values(*column))})
        response['Access-Control-Allow-Origin'] = '*'
        return response


SOURCE_RELATED_FIELDS = {
    'Variant_in_ENIGMA': 'enigma_reports__isnull',
    'Variant_in_ClinVar': 'clinvar_data__isnull',
    'Variant_in_LOVD': 'lovd_data__isnull',
    'Variant_in_GnomAD': 'gnomad_data__isnull',
}

ALL_SOURCES = list(SOURCE_RELATED_FIELDS.keys())


def apply_sources(query, include, exclude):
    # if there are multiple sources given then OR them:
    # the row must match in at least one column
    if len(include) > 0:
        if include == ['all']:
            include = ALL_SOURCES
        include_list = (Q(**{SOURCE_RELATED_FIELDS[s]: False}) for s in include if s in SOURCE_RELATED_FIELDS)
        query = query.filter(reduce(__or__, include_list)).distinct()
    else:
        # exclude all sources if none are included
        for source in (exclude or ALL_SOURCES):
            if source in SOURCE_RELATED_FIELDS:
                query = query.filter(**{SOURCE_RELATED_FIELDS[source]: True})
    return query


def normalize_filter_values(filterValues):
    return [fV.replace('Likely Benign', 'Likely benign').replace('Likely Pathogenic', 'Likely pathogenic') for fV in filterValues]


def apply_filters(query, filterValues, filters, quotes=''):
    # if there are multiple filters the row must match all the filters
    normalizedFilterValues = normalize_filter_values(filterValues)
    for column, value in zip(filters, normalizedFilterValues):
        if column == 'id':
            query = query.filter(**{column: value})
        else:
            query = query.extra(
                where=["\"{0}\" LIKE %s".format(column)],
                params=["{0}{1}%{0}".format(quotes, value)]
            )
    return query


def add_paren_to_hgvs_protein_if_absent(value):
    if value.startswith('p.') and '(' not in value:
        return value[:2] + '(' + value[2:]
    else:
        return value


def remove_disallowed_chars(search_term):
    for disallowed in DISALLOWED_SEARCH_CHARS:
        search_term = search_term.replace(disallowed, '')
    return search_term


def apply_search(query, search_term, quotes=''):
    '''
    NOTE: there is some additional handling of search terms on the front-end in
    website/js/hgvs.js. hgvs.js methods are called before sending the query to the
    backend.

    Below are examples of all special case searches that don't match our data schema
    but are handled in this method. Each example contains a user submitted search
    followed by the colon delimited fields that they represent (note that some fields
    contain colons in and of themselves, e.g. Genomic Coordinates).

        User submitted search --> Field:Field

        BRCA1:958C>G --> Gene_Symbol:BIC_Nomenclature
        BRCA1:c.839C>G --> Gene_Symbol:HGVS_cDNA
        NM_007294.3:958C>G --> Reference_Sequence:BIC_Nomenclature
        NM_007294.3:c.839C>G --> Reference_Sequence:HGVS_cDNA
        BRCA1:p.(Ala280Gly) --> Gene_Symbol:HGVS_Protein.split(':')[1] (HGVS_Protein is actually stored as NP_009225.1:p.(Ala280Gly), so this has to be split on the ":")
        BRCA1:A280G --> Gene_Symbol:Protein_Change
        NP_009225.1:p.(Ala280Gly) --> HGVS_Protein
        NP_009225.1:A280G --> HGVS_Protein.split(':')[0]:Protein_Change
    '''
    search_term = search_term.lower().strip()
    search_term = remove_disallowed_chars(search_term)
    clinvar_accession = False

    # Accept search by ClinGen Allele Registry id (CA_ID)
    if search_term.startswith('ca'):
        return query.filter(Q(CA_ID__icontains=search_term)), 0

    # Accept search by ga4gh VR id
    if search_term.startswith('ga4gh'):
        return query.filter(Q(VRS__digest__icontains=search_term)), 0

    # Accept only full clinvar accession numbers
    if search_term.startswith('scv') and len(search_term) >= 12:
        clinvar_accession = True

    p_hgvs_protein_colon = re.compile("^np_[0-9]{6}.[0-9]:")
    m_hgvs_protein_colon = p_hgvs_protein_colon.match(search_term)
    p_hgvs_protein_space = re.compile("^np_[0-9]{6}.[0-9] ")
    m_hgvs_protein_space = p_hgvs_protein_space.match(search_term)

    p_reference_sequence_colon = re.compile("^nm_[0-9]{6}.[0-9]:")
    m_reference_sequence_colon = p_reference_sequence_colon.match(search_term)
    p_reference_sequence_space = re.compile("^nm_[0-9]{6}.[0-9] ")
    m_reference_sequence_space = p_reference_sequence_space.match(search_term)

    has_gene_symbol_prefix = False
    for accepted_prefix in ['brca1:', 'brca2:', 'brca1 ', 'brca2 ']:
        if search_term.startswith(accepted_prefix):
            has_gene_symbol_prefix = True

    # Handle HGVS_Protein searches
    if m_hgvs_protein_space or m_hgvs_protein_colon:
        prefix = search_term[:11]
        suffix = search_term[12:]

        # accept hgvs_protein sequences without parentheses
        suffix = add_paren_to_hgvs_protein_if_absent(suffix)

        # values in synonyms column are separated by commas
        comma_prefixed_suffix = ',' + suffix
        results = query.filter(HGVS_Protein__istartswith=prefix).filter(
            Q(Protein_Change__istartswith=suffix) |
            Q(HGVS_Protein__icontains=suffix) |
            Q(Synonyms__icontains=comma_prefixed_suffix) |
            Q(Synonyms__istartswith=suffix)
        ) | query.filter(Q(HGVS_Protein__icontains=search_term) | Q(Synonyms__icontains=search_term))
        non_synonyms = results.filter(Protein_Change__istartswith=suffix) | query.filter(HGVS_Protein__icontains=search_term)

    # Handle gene symbol prefixed searches
    elif has_gene_symbol_prefix:
        prefix = search_term[:5]
        suffix = search_term[6:]

        suffix = add_paren_to_hgvs_protein_if_absent(suffix)

        comma_prefixed_suffix = ',' + suffix
        # need to check synonym column for colon prefixes in the case of HGVS_cDNA and HGVS_Protein fields
        colon_prefixed_suffix = ':' + suffix
        results = query.filter(Gene_Symbol__iexact=prefix).filter(
            Q(HGVS_cDNA__icontains=suffix) |
            Q(HGVS_Protein__icontains=suffix) |
            Q(BIC_Nomenclature__istartswith=suffix) |
            Q(Protein_Change__istartswith=suffix) |
            Q(Synonyms__icontains=comma_prefixed_suffix) |
            Q(Synonyms__icontains=colon_prefixed_suffix) |
            Q(Synonyms__istartswith=suffix)
        ) | query.filter(Synonyms__icontains=search_term)
        non_synonyms = results.filter(
            Q(HGVS_cDNA__icontains=suffix) |
            Q(HGVS_Protein__icontains=suffix) |
            Q(BIC_Nomenclature__istartswith=suffix) |
            Q(Protein_Change__istartswith=suffix)
        )

    # Handle Reference_Sequence prefixed searches
    elif m_reference_sequence_space or m_reference_sequence_colon:
        prefix = search_term[:11]
        suffix = search_term[12:]
        comma_prefixed_suffix = ',' + suffix
        results = query.filter(Reference_Sequence__iexact=prefix).filter(
            Q(HGVS_cDNA__icontains=suffix) |
            Q(BIC_Nomenclature__istartswith=suffix) |
            Q(Synonyms__icontains=comma_prefixed_suffix) |
            Q(Synonyms__istartswith=suffix)
        ) | query.filter(Synonyms__icontains=search_term)
        non_synonyms = results.filter(
            Q(HGVS_cDNA__icontains=suffix) |
            Q(BIC_Nomenclature__istartswith=suffix)
        )

    # Handle clinvar accession numbers
    elif clinvar_accession is True:
        results = query.filter(
            Q(clinvar_data__clinvar_reports__SCV__icontains=search_term) |
            Q(enigma_reports__ClinVarAccession__icontains=search_term)
        ).distinct()
        non_synonyms = results

    # Generic searches (no prefixes)
    else:
        results = query.filter(
            Q(Pathogenicity__icontains=search_term) |
            Q(Synonyms__icontains=search_term) |
            Q(Gene_Symbol__icontains=search_term) |
            Q(HGVS_cDNA__icontains=search_term) |
            Q(BIC_Nomenclature__icontains=search_term) |
            Q(HGVS_Protein__icontains=search_term) |
            Q(Protein_Change__icontains=search_term)
        )

        non_synonyms = query.filter(
            Q(Pathogenicity__icontains=search_term) |
            Q(Gene_Symbol__icontains=search_term) |
            Q(HGVS_cDNA__icontains=search_term) |
            Q(BIC_Nomenclature__icontains=search_term) |
            Q(HGVS_Protein__icontains=search_term) |
            Q(Protein_Change__icontains=search_term)
        )

    synonyms_count = results.count() - non_synonyms.count()

    return results, synonyms_count


def apply_order(query, order_by, direction):
    # special case for HGVS columns
    if order_by in ('HGVS_cDNA', 'HGVS_Protein'):
        order_by = 'Gene_Symbol'
    if direction == 'descending':
        order_by = '-' + order_by
    return query.order_by(order_by, 'Pathogenicity')


def select_page(query, page_size, page_num):
    if page_size:
        start = page_size * page_num
        end = start + page_size
        return query[start:end]
    return query


def autocomplete(request):
    cursor = connection.cursor()
    term = request.GET.get('term')

    '''If a release is specified in the query, only return autocomplete
    suggestions for specified release, otherwise default to suggestions
    for the latest release'''
    if 'release' in request.GET:
        release = request.GET.get('release')
    else:
        cursor.execute("""SELECT MAX(id) FROM data_release""")
        release = cursor.fetchone()[0]

    limit = int(request.GET.get('limit', 10))

    cursor.execute(
        """SELECT word FROM words
        WHERE word LIKE %s
        AND char_length(word) >= 3
        AND release_id = %s
        ORDER BY word""",
        ["%s%%" % term, release])

    rows = cursor.fetchall()

    response = JsonResponse({'suggestions': rows[:limit]})
    response['Access-Control-Allow-Origin'] = '*'
    return response


def sanitise_term(term):
    # Escape all non alphanumeric characters
    term = re.escape(term)
    # Enable prefix search
    term += ":*"
    return term
