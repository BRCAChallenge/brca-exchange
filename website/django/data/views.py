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
    Variant, VariantDiff, CurrentVariant, DataRelease, ChangeType, Report, ReportDiff,
    InSilicoPriors, VariantPaper, Paper, VariantRepresentation
)
from django.views.decorators.http import require_http_methods

import google.protobuf.json_format as json_format
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
    change_types = {x['name']: x['id'] for x in list(ChangeType.objects.values())}
    for release in releases:
        variants = Variant.objects.filter(Data_Release_id=release['id'])
        release['variants_added'] = variants.filter(Change_Type_id=change_types['new']).count()
        release['variants_classified'] = variants.filter(
            Q(Change_Type_id=change_types['changed_classification']) | Q(Change_Type_id=change_types['added_classification'])).count()
        release['variants_modified'] = variants.filter(
            Q(Change_Type_id=change_types['added_information']) | Q(Change_Type_id=change_types['changed_information'])).count()
        release['variants_deleted'] = variants.filter(Change_Type_id=change_types['deleted']).count()
    response = JsonResponse({"releases": list(releases), "latest": latest})
    response['Access-Control-Allow-Origin'] = '*'
    return response

def variant_counts(request):
    query = CurrentVariant.objects.all().exclude(Change_Type__name='deleted')
    total_count = query.count()
    brca1_count = query.filter(Gene_Symbol='BRCA1').count()
    brca2_count = query.filter(Gene_Symbol='BRCA2').count()
    query = query.filter(Variant_in_ENIGMA=True)
    enigma_count = query.count()
    enigma_pathogenic_count = query.filter(Pathogenicity_expert='Pathogenic').count()
    enigma_benign_count = query.filter(Pathogenicity_expert__contains='Benign').count()
    enigma_likely_benign_count = query.filter(Pathogenicity_expert__contains='Likely benign').count()
    enigma_likely_pathogenic_count = query.filter(Pathogenicity_expert__contains='Likely pathogenic').count()
    query_brca1 = query.filter(Gene_Symbol='BRCA1')
    brca1_enigma_pathogenic_count = query_brca1.filter(Pathogenicity_expert='Pathogenic').count()
    brca1_enigma_benign_count = query_brca1.filter(Pathogenicity_expert__contains='Benign').count()
    brca1_enigma_likely_benign_count = query_brca1.filter(Pathogenicity_expert__contains='Likely benign').count()
    brca1_enigma_likely_pathogenic_count = query_brca1.filter(Pathogenicity_expert__contains='Likely pathogenic').count()
    query_brca2 = query.filter(Gene_Symbol='BRCA2')
    brca2_enigma_pathogenic_count = query_brca2.filter(Pathogenicity_expert='Pathogenic').count()
    brca2_enigma_benign_count = query_brca2.filter(Pathogenicity_expert__contains='Benign').count()
    brca2_enigma_likely_benign_count = query_brca2.filter(Pathogenicity_expert__contains='Likely benign').count()
    brca2_enigma_likely_pathogenic_count = query_brca2.filter(Pathogenicity_expert__contains='Likely pathogenic').count()
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
    variant_id = int(request.GET.get('variant_id'))

    variant = Variant.objects.get(id=variant_id)
    key = variant.Genomic_Coordinate_hg38
    query = Variant.objects.filter(Genomic_Coordinate_hg38=key)\
        .order_by('-Data_Release_id')\
        .select_related('Data_Release')\
        .select_related('Mupit_Structure')\
        .select_related('insilicopriors')

    variant_versions = list(map(variant_to_dict, query))
    response = JsonResponse({"data": variant_versions})
    response['Access-Control-Allow-Origin'] = '*'
    return response


def variantreps(request):
    vr_reps = list(
        VariantRepresentation.objects.raw("""
        select VR.id, VR."Genomic_Coordinate_hg38", CV.id as Variant_id, VR."Description" from data_variantrepresentation VR
        inner join currentvariant CV on CV."Genomic_Coordinate_hg38" = VR."Genomic_Coordinate_hg38"
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
    variants = CurrentVariant.objects.values_list('id')
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


def variant_reports(request, variant_id):
    variant_id = int(variant_id)
    query = Report.objects.filter(Variant_id=variant_id)
    report_versions = []
    for report in query:
        key = None
        if report.Source == "ClinVar":
            key = report.SCV_ClinVar
            if not key or key == '-':
                # if no key is available, skip report history
                report_query = [report]
            else:
                # extend the selection w/reports that have matching keys,
                # but only up until the requested variants' release
                report_query = Report.objects\
                    .filter(Data_Release_id__lte=report.Data_Release.id, SCV_ClinVar=key)\
                    .order_by('-Data_Release_id').select_related('Data_Release')
            report_versions.extend(list(map(report_to_dict, report_query)))
        elif report.Source == "LOVD":
            key = report.Submission_ID_LOVD
            if not key or key == '-':
                # if no key is available, skip report history
                report_query = [report]
            else:
                # extend the selection w/reports that have matching keys,
                # but only up until the requested variants' release (i.e., same as for ClinVar)
                report_query = Report.objects\
                    .filter(Data_Release_id__lte=report.Data_Release.id, Submission_ID_LOVD=key)\
                    .order_by('-Data_Release_id').select_related('Data_Release')
            report_versions.extend(list(map(report_to_dict, report_query)))

    response = JsonResponse({"data": report_versions})
    response['Access-Control-Allow-Origin'] = '*'
    return response

def variant_papers(request):
    variant_id = int(request.GET.get('variant_id'))
    variant = Variant.objects.get(id=variant_id)
    variant_name = variant.Genomic_Coordinate_hg38
    variantpapers = VariantPaper.objects.select_related('paper').filter(variant_hg38=variant_name).all()
    for variantpaper in variantpapers:
        # year of 0000 means year could not be found during a crawl
        if variantpaper.paper.year == 0000:
            variantpaper.paper.year = "Unknown"
    variantpapers = [dict(model_to_dict(vp.paper), **{"mentions": vp.mentions, "points": vp.points}) for vp in variantpapers]
    response = JsonResponse({"data": variantpapers}, safe=False)
    response['Access-Control-Allow-Origin'] = '*'
    return response

def variant_to_dict(variant_object):
    change_types_map = {x['name']:x['id'] for x in list(ChangeType.objects.values())}
    variant_dict = model_to_dict(variant_object)
    variant_dict["Data_Release"] = model_to_dict(variant_object.Data_Release)
    if variant_object.Mupit_Structure is not None:
        variant_dict["Mupit_Structure"] = model_to_dict(variant_object.Mupit_Structure)
    variant_dict["Data_Release"]["date"] = variant_object.Data_Release.date
    variant_dict["Change_Type"] = ChangeType.objects.get(id=variant_dict["Change_Type"]).name

    try:
        variant_dict["priors"] = model_to_dict(variant_object.insilicopriors)
    except InSilicoPriors.DoesNotExist:
        variant_dict["priors"] = None

    try:
        variant_diff = VariantDiff.objects.get(variant_id=variant_object.id)
        variant_dict["Diff"] = variant_diff.diff
    except VariantDiff.DoesNotExist:
        variant_dict["Diff"] = None
    return variant_dict


def report_to_dict(report_object):
    change_types_map = {x['name']:x['id'] for x in list(ChangeType.objects.values())}
    report_dict = model_to_dict(report_object)
    report_dict["Data_Release"] = model_to_dict(report_object.Data_Release)
    report_dict["Data_Release"]["date"] = report_object.Data_Release.date
    report_dict["Change_Type"] = ChangeType.objects.get(id=report_dict["Change_Type"]).name

    if report_object.Source == "ClinVar":
        # don't display ClinVar report diffs prior to April 2018
        cutoff_date = datetime.strptime('Apr 1 2018  12:00AM', '%b %d %Y %I:%M%p')
    elif report_object.Source == "LOVD":
        # don't display LOVD report diffs prior to November 4 2018 (we
        # updated the definition of LOVD submissions in the early November
        # release, so it only makes sense to show diffs from following releases)
        cutoff_date = datetime.strptime('Nov 4 2018  12:00AM', '%b %d %Y %I:%M%p')
    try:
        if report_dict["Data_Release"]["date"] < cutoff_date:
            report_dict["Diff"] = None
            return report_dict
    except UnboundLocalError as e:
        logging.error(repr(e))

    try:
        report_diff = ReportDiff.objects.get(report_id=report_object.id)
        report_dict["Diff"] = report_diff.report_diff
    except ReportDiff.DoesNotExist:
        report_dict["Diff"] = None


    return report_dict


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
    release = request.GET.get('release')
    change_types = request.GET.getlist('change_types')
    change_types_map = {x['name']:x['id'] for x in list(ChangeType.objects.values())}
    show_deleted = (request.GET.get('show_deleted', False) != False)
    deleted_count = 0
    synonyms_count = 0
    release_name = None

    if release:
        query = Variant.objects.filter(Data_Release_id=int(release))
        if(change_types):
            change_types = [change_types_map[c] for c in [c for c in change_types if c in change_types_map]]
            query = query.filter(Change_Type_id__in=change_types)
        release_name = DataRelease.objects.get(id=int(release)).name
    else:
        query = CurrentVariant.objects

    if format == 'csv' or format == 'tsv':
        quotes = '\''
    else:
        quotes = ''

    query = apply_sources(query, include, exclude)

    if filters:
        query = apply_filters(query, filter_values, filters, quotes=quotes)

    if search_term:
        query, synonyms_count = apply_search(query, search_term, quotes=quotes, release=release)

    if not show_deleted and not release:
        deleted_count = query.filter(Change_Type_id=change_types_map['deleted']).count()
        query = query.exclude(Change_Type_id=change_types_map['deleted'])

    if order_by:
        query = apply_order(query, order_by, direction)

    if format == 'csv' or format == 'tsv':
        cursor = connection.cursor()

        # create an in-memory StringIO object that we can use to buffer the results from the server
        # (it'll get released when it goes out of scope, so unlike a real file we don't need to close it)
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
        # call list() now to evaluate the query
        response = JsonResponse({'count': count, 'deletedCount': deleted_count, 'synonyms': synonyms_count, 'releaseName': release_name, 'data': list(query.values(*column))})
        response['Access-Control-Allow-Origin'] = '*'
        return response


def apply_sources(query, include, exclude):
    # if there are multiple sources given then OR them:
    # the row must match in at least one column
    if len(include) > 0:
        include_list = (Q(**{column: True}) for column in include)
        query = query.filter(reduce(__or__, include_list))
    else:
        # exclude all sources if none are included
        exclude = [f.name for f in Variant._meta.get_fields() if "Variant_in" in f.name]
    if exclude:
        exclude_dict = {exclusion: False for exclusion in exclude}
        query = query.filter(**exclude_dict)
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


def apply_search(query, search_term, quotes='', release=None):
    '''
    NOTE: there is some additional handling of search terms on the front-end in
    website/js/hgvs.js. hgvs.js methods are called before sending the query to the
    backend.

    Below are examples of all special case searches that don't match our data schema
    but are handled in this method. Each example contains a user submitted search
    followed by the colon delimited fields that they represent (note that some fields
    contain colons in and of themselves, e.g. Genomic Coordinates).

        User submitted search --> Field:Field

        BRCA1:chr17:g.43094692:G>C --> Gene_Symbol:Genomic_Coordinate_hg38
        BRCA1:chr17:g.41246709:G>C --> Gene_Symbol:Genomic_Coordinate_hg37
        BRCA1:chr17:g.38500235:G>C --> Gene_Symbol:Genomic_Coordinate_hg36
        BRCA1:958C>G --> Gene_Symbol:BIC_Nomenclature
        BRCA1:c.839C>G --> Gene_Symbol:HGVS_cDNA
        NM_007294.3:chr17:g.43094692:G>C --> Reference_Sequence:Genomic_Coordinate_hg38
        NM_007294.3:chr17:g.41246709:G>C --> Reference_Sequence:Genomic_Coordinate_hg37
        NM_007294.3:chr17:g.38500235:G>C --> Reference_Sequence:Genomic_Coordinate_hg36
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

    # Accept only full clinvar accession numbers
    if search_term.startswith('scv') and len(search_term) >= 12:
        clinvar_accession = True

    # Accept genomic coordinates with or without a 'g.' before the position
    if 'chr17:' in search_term and 'g.' not in search_term:
        search_term = search_term.replace('chr17:', 'chr17:g.')
    if 'chr13:' in search_term and 'g.' not in search_term:
        search_term = search_term.replace('chr13:', 'chr13:g.')

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
            Q(Genomic_Coordinate_hg38__istartswith=suffix) |
            Q(Genomic_Coordinate_hg37__istartswith=suffix) |
            Q(Genomic_Coordinate_hg36__istartswith=suffix) |
            Q(BIC_Nomenclature__istartswith=suffix) |
            Q(Protein_Change__istartswith=suffix) |
            Q(Synonyms__icontains=comma_prefixed_suffix) |
            Q(Synonyms__icontains=colon_prefixed_suffix) |
            Q(Synonyms__istartswith=suffix)
        ) | query.filter(Synonyms__icontains=search_term)
        non_synonyms = results.filter(
            Q(HGVS_cDNA__icontains=suffix) |
            Q(HGVS_Protein__icontains=suffix) |
            Q(Genomic_Coordinate_hg38__istartswith=suffix) |
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
            Q(Genomic_Coordinate_hg38__istartswith=suffix) |
            Q(Genomic_Coordinate_hg37__istartswith=suffix) |
            Q(Genomic_Coordinate_hg36__istartswith=suffix) |
            Q(BIC_Nomenclature__istartswith=suffix) |
            Q(Synonyms__icontains=comma_prefixed_suffix) |
            Q(Synonyms__istartswith=suffix)
        ) | query.filter(Synonyms__icontains=search_term)
        non_synonyms = results.filter(
            Q(HGVS_cDNA__icontains=suffix) |
            Q(Genomic_Coordinate_hg38__istartswith=suffix) |
            Q(BIC_Nomenclature__istartswith=suffix)
        )
    # Handle clinvar accession numbers
    elif clinvar_accession is True:
        results = query.filter(
            Q(SCV_ClinVar__icontains=search_term) |
            Q(ClinVarAccession_ENIGMA__icontains=search_term)
        )

        # filter against synonym fields
        non_synonyms = query.filter(
            Q(SCV_ClinVar__icontains=search_term) |
            Q(ClinVarAccession_ENIGMA__icontains=search_term)
        )

        # Generic searches (no prefixes)
    else:
        # filter non-special-case searches against the following fields
        results = query.filter(
            Q(Pathogenicity_expert__icontains=search_term) |
            Q(Genomic_Coordinate_hg38__icontains=search_term) |
            Q(Genomic_Coordinate_hg37__icontains=search_term) |
            Q(Genomic_Coordinate_hg36__icontains=search_term) |
            Q(Synonyms__icontains=search_term) |
            Q(Gene_Symbol__icontains=search_term) |
            Q(HGVS_cDNA__icontains=search_term) |
            Q(BIC_Nomenclature__icontains=search_term) |
            Q(HGVS_Protein__icontains=search_term) |
            Q(Protein_Change__icontains=search_term)
        )

        # filter against synonym fields
        non_synonyms = query.filter(
            Q(Pathogenicity_expert__icontains=search_term) |
            Q(Genomic_Coordinate_hg38__icontains=search_term) |
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
        order_by = 'Genomic_Coordinate_hg38'
    if direction == 'descending':
        order_by = '-' + order_by
    return query.order_by(order_by, 'Pathogenicity_expert')


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
