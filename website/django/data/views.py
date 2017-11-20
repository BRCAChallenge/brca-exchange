import os
import re
import tempfile
import json
from operator import __or__
from django.db import connection
from django.db.models import Q
from django.db.models import Value
from django.db.models.functions import Concat
from django.forms.models import model_to_dict
from django.http import JsonResponse, HttpResponse, HttpResponseBadRequest
from django.views.decorators.gzip import gzip_page

from .models import Variant, VariantDiff, CurrentVariant, DataRelease, ChangeType
from django.views.decorators.http import require_http_methods

# GA4GH related imports
from ga4gh.schemas.ga4gh import variant_service_pb2 as variant_service
from ga4gh.schemas.ga4gh import variants_pb2 as variants
from ga4gh.schemas.ga4gh import metadata_service_pb2 as metadata_service
from ga4gh.schemas.ga4gh import metadata_pb2 as metadata

import google.protobuf.json_format as json_format


def releases(request):
    release_id = request.GET.get('release_id')
    if release_id:
        releases = DataRelease.objects.filter(id=release_id).values().all()
    else:
        releases = DataRelease.objects.values().all()
    latest = DataRelease.objects.order_by('-id')[0].id
    change_types = {x['name']: x['id'] for x in ChangeType.objects.values().all()}
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
    query_brca1 = query.filter(Gene_Symbol='BRCA1')
    brca1_enigma_pathogenic_count = query_brca1.filter(Pathogenicity_expert='Pathogenic').count()
    brca1_enigma_benign_count = query_brca1.filter(Pathogenicity_expert__contains='Benign').count()
    brca1_enigma_likely_benign_count = query_brca1.filter(Pathogenicity_expert__contains='Likely benign').count()
    query_brca2 = query.filter(Gene_Symbol='BRCA2')
    brca2_enigma_pathogenic_count = query_brca2.filter(Pathogenicity_expert='Pathogenic').count()
    brca2_enigma_benign_count = query_brca2.filter(Pathogenicity_expert__contains='Benign').count()
    brca2_enigma_likely_benign_count = query_brca2.filter(Pathogenicity_expert__contains='Likely benign').count()
    response = JsonResponse({
        "total": total_count,
        "brca1": {
            "total": brca1_count,
            "pathogenic": brca1_enigma_pathogenic_count,
            "benign": brca1_enigma_benign_count,
            "likelyBenign": brca1_enigma_likely_benign_count },
        "brca2": {
            "total": brca2_count,
            "pathogenic": brca2_enigma_pathogenic_count,
            "benign": brca2_enigma_benign_count,
            "likelyBenign": brca2_enigma_likely_benign_count },
        "enigma": enigma_count,
        "enigmaPathogenic": enigma_pathogenic_count,
        "enigmaBenign": enigma_benign_count,
        "enigmaLikelyBenign": enigma_likely_benign_count })
    response['Access-Control-Allow-Origin'] = '*'
    return response

def variant(request):
    variant_id = int(request.GET.get('variant_id'))

    variant = Variant.objects.get(id=variant_id)
    key = variant.Genomic_Coordinate_hg38

    query = Variant.objects.filter(Genomic_Coordinate_hg38=key).order_by('-Data_Release_id').select_related('Data_Release')

    variant_versions = map(variant_to_dict, query)
    response = JsonResponse({"data": variant_versions})
    response['Access-Control-Allow-Origin'] = '*'
    return response


def variant_to_dict(variant_object):
    change_types_map = {x['name']:x['id'] for x in ChangeType.objects.values().all()}
    variant_dict = model_to_dict(variant_object)
    variant_dict["Data_Release"] = model_to_dict(variant_object.Data_Release)
    variant_dict["Data_Release"]["date"] = variant_object.Data_Release.date
    variant_dict["Change_Type"] = ChangeType.objects.get(id=variant_dict["Change_Type"]).name
    try:
        variant_diff = VariantDiff.objects.get(variant_id=variant_object.id)
        variant_dict["Diff"] = variant_diff.diff
    except VariantDiff.DoesNotExist:
        print "Variant Diff does not exist for Variant", variant_object.Genomic_Coordinate_hg38, "from release", variant_object.Data_Release.id
        variant_dict["Diff"] = None
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
    release = request.GET.get('release')
    change_types = request.GET.getlist('change_types')
    change_types_map = {x['name']:x['id'] for x in ChangeType.objects.values().all()}
    show_deleted = (request.GET.get('show_deleted', False) != False)
    deleted_count = 0
    synonyms_count = 0

    if release:
        query = Variant.objects.filter(Data_Release_id=int(release))
        if(change_types):
            change_types = map(lambda c: change_types_map[c], filter(lambda c: c in change_types_map, change_types))
            query = query.filter(Change_Type_id__in=change_types)
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

    if format == 'csv':

        cursor = connection.cursor()
        with tempfile.NamedTemporaryFile() as f:
            os.chmod(f.name, 0606)
            query = "COPY ({}) TO '{}' WITH DELIMITER ',' CSV HEADER".format(query.query, f.name)
            # HACK to add quotes around search terms
            query = re.sub(r'LIKE UPPER\((.+?)\)', r"LIKE UPPER('\1')", query)
            cursor.execute(query)

            response = HttpResponse(f.read(), content_type='text/csv')
            response['Content-Disposition'] = 'attachment;filename="variants.csv"'
            return response

    elif format == 'tsv':

        cursor = connection.cursor()
        with tempfile.NamedTemporaryFile() as f:
            os.chmod(f.name, 0606)
            query = "COPY ({}) TO '{}' WITH DELIMITER '\t' CSV HEADER".format(query.query, f.name)
            # HACK to add quotes around search terms
            query = re.sub(r'LIKE UPPER\((.+?)\)', r"LIKE UPPER('\1')", query)
            cursor.execute(query)

            response = HttpResponse(f.read(), content_type='text/csv')
            response['Content-Disposition'] = 'attachment;filename="variants.tsv"'
            return response

    elif format == 'json':
        count = query.count()
        query = select_page(query, page_size, page_num)
        # call list() now to evaluate the query
        response = JsonResponse({'count': count, 'deletedCount': deleted_count, 'synonyms': synonyms_count, 'data': list(query.values(*column))})
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


def apply_filters(query, filterValues, filters, quotes=''):
    # if there are multiple filters the row must match all the filters
    for column, value in zip(filters, filterValues):
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


@require_http_methods(["POST"])
def search_variants(request):
    """Handles requests to the /variants/search method"""
    conditional = validate_search_variants_request(request)
    if conditional:
        return conditional
    else:
        try:
            protocol_variable = json_format.Parse(request.body, variant_service.SearchVariantsRequest())
        except Exception as e:
            return HttpResponseBadRequest(json.dumps({"message": e.message.replace("\"", "'")}),
                                          content_type='application/json')
        variant_set_id = protocol_variable.variant_set_id
        reference_name = protocol_variable.reference_name
        start = protocol_variable.start
        end = protocol_variable.end
        page_size = protocol_variable.page_size
        page_token = protocol_variable.page_token
    if not page_size or page_size == 0:
        page_size = DEFAULT_PAGE_SIZE
    if not page_token:
        page_token = '0'

    response = variant_service.SearchVariantsResponse()
    variants = CurrentVariant.objects
    dataset_id, reference_genome = variant_set_id.split('-')
    if dataset_id != DATASET_ID:
        return HttpResponseBadRequest(
                json.dumps(ErrorMessages['variantSetId']),
                content_type='application/json')
    variants = range_filter(reference_genome, variants, reference_name, start, end)
    variants = ga4gh_brca_page(variants, int(page_size), int(page_token))

    ga_variants = []
    for i in variants.values():
        try:
            ga_variants.append(brca_to_ga4gh(i, reference_genome))
        except ValueError as e:
            print e
    if len(ga_variants) > page_size:
        ga_variants.pop()
        page_token = str(1 + int(page_token))
        response.next_page_token = page_token

    response.variants.extend(ga_variants)
    resp = json_format.MessageToDict(response, True)
    return JsonResponse(resp)


def range_filter(reference_genome, variants, reference_name, start, end):
    """Filters variants by range depending on the reference_genome"""
    if 'chr' in reference_name:
        reference_name = reference_name.replace('chr', '')
    variants = variants.filter(Chr=reference_name)
    if reference_genome == 'hg36':
        variants = variants.order_by('Hg36_Start')
        variants = variants.filter(Hg36_Start__lt=end, Hg36_End__gt=start)
    elif reference_genome == 'hg37':
        variants = variants.order_by('Hg37_Start')
        variants = variants.filter(Hg37_Start__lt=end, Hg37_End__gt=start)
    elif reference_genome == 'hg38':
        variants = variants.order_by('Hg38_Start')
        variants = variants.filter(Hg38_Start__lt=end, Hg38_End__gt=start)
    return variants

def ga4gh_brca_page(query, page_size, page_token):
    """Filters django queries by page for GA4GH requests"""
    start = page_size * page_token
    end = start + page_size + 1
    return query[start:end]

def brca_to_ga4gh(brca_variant, reference_genome):
    """Function that translates elements in BRCA-database to GA4GH format."""
    brca_variant = {k: unicode(v).encode("utf-8") for k,v in brca_variant.iteritems()}
    variant = variants.Variant()
    bases = brca_variant['Genomic_Coordinate_' + reference_genome].split(':')[2]
    variant.reference_bases, alternbases = bases.split('>')
    for i in alternbases.split(","):
        variant.alternate_bases.append(i)
    variant.created = 0
    variant.updated = 0
    variant.reference_name = brca_variant['Chr']
    if reference_genome == 'hg36':
        variant.start = int(brca_variant['Hg36_Start'])
        variant.end = int(brca_variant['Hg36_End'])
    elif reference_genome == 'hg37':
        variant.start = int(brca_variant['Hg37_Start'])
        variant.end = int(brca_variant['Hg37_End'])
    elif reference_genome == 'hg38':
        variant.start = int(brca_variant['Hg38_Start'])
        variant.end = int(brca_variant['Hg38_End'])
    variant.id = '{}-{}'.format(reference_genome, str(brca_variant['id']))
    variant.variant_set_id = '{}-{}'.format(DATASET_ID, reference_genome)
    names = brca_variant['Synonyms'].split(',')
    for name in names:
        variant.names.append(name)
    for key in brca_variant:
        if brca_variant[key] != '-' and brca_variant[key] != '':
            variant.info[str(key)].append(brca_variant[key])
    return variant

def validate_search_variants_request(request):
    """Auxiliary function which validates search variants requests"""
    if not request.body:
        return HttpResponseBadRequest(
            json.dumps(ErrorMessages['emptyBody']),
            content_type='application/json')
    else:
        request_dict = json.loads(request.body)
        if not request_dict.get('variantSetId'):
            return HttpResponseBadRequest(
                json.dumps(ErrorMessages['variantSetId']),
                content_type='application/json')
        elif not request_dict.get('referenceName'):
            return HttpResponseBadRequest(
                json.dumps(ErrorMessages['referenceName']),
                content_type='application/json')
        elif not request_dict.get('start'):
            return HttpResponseBadRequest(
                json.dumps(ErrorMessages['start']),
                content_type='application/json')
        elif not request_dict.get('end'):
            return HttpResponseBadRequest(
                json.dumps(ErrorMessages['end']),
                content_type='application/json')
        elif int(request_dict.get('start')) >= int(request_dict.get('end')):
            return HttpResponseBadRequest(
                json.dumps(ErrorMessages['invalidPositions']),
                content_type='application/json')
        else:
                # Make sure the variant set ID is well formed
            ids = request_dict.get('variantSetId').split('-')
            reference_name = request_dict.get('referenceName')
            # A malformed ID
            if len(ids) < 2:
                return HttpResponse(
                    json.dumps(ErrorMessages['variantSetId']),
                               content_type='application/json',
                               status=404)
            reference_genome = ids[1]
            if reference_genome not in SET_IDS:
                return HttpResponse(
                    json.dumps(ErrorMessages['variantSetId']),
                    content_type='application/json',
                    status=404)
            return None

def validate_search_variant_sets_request(request):
    """Auxiliary function which validates search variant sets requests"""
    if not request.body:
        return HttpResponseBadRequest(
            json.dumps(ErrorMessages['emptyBody']),
            content_type='application/json')
    else:
        request_dict = json.loads(request.body)
        if not request_dict.get('datasetId'):
            return HttpResponseBadRequest(
                json.dumps(ErrorMessages['datasetId']),
                content_type='application/json')
        else:
            return None

@require_http_methods(['GET'])
def get_variant(request, variant_id):
    """Handles requests to the /variants/<variant id> endpoint"""
    if not variant_id:
        return HttpResponseBadRequest(
            json.dumps(ErrorMessages['variantId']),
            content_type='application/json')
    else:
        set_id, v_id = variant_id.split('-')
        if set_id in SET_IDS:
            variants = CurrentVariant.objects.values()
            try:
                variant = variants.get(id=int(v_id))
            except Exception:
                return HttpResponse(
                    json.dumps(ErrorMessages['notFoundId']),
                    content_type='application/json',
                    status=404)
            ga_variant = brca_to_ga4gh(variant, set_id)
            response = json_format.MessageToDict(ga_variant, True)
            return JsonResponse(response)
        else:
            return HttpResponse(
                json.dumps(ErrorMessages['notFoundId']),
                content_type='application/json',
                status=404)

@require_http_methods(['POST'])
def search_variant_sets(request):
    """Handles requests at the /variantsets/search endpoint"""
    invalid_request = validate_search_variant_sets_request(request)
    if invalid_request:
        return invalid_request
    else:
        try:
            req_dict = json_format.Parse(request.body, variant_service.SearchVariantSetsRequest())
        except Exception as e :
            return HttpResponseBadRequest(json.dumps({"message": e.message.replace("\"", "'")}),
                                          content_type='application/json')
        dataset_id = req_dict.dataset_id
        page_size = req_dict.page_size
        page_token = req_dict.page_token
        if dataset_id != DATASET_ID:
            """Bad Request returns empty response"""
            return JsonResponse(
                json_format.MessageToDict(
                    variant_service.SearchCallSetsResponse(), True))
    if not page_size or page_size == 0:
        page_size = DEFAULT_PAGE_SIZE
    if not page_token:
        page_token = '0'

    response = variant_service.SearchVariantSetsResponse()
    variant_sets_list = [obtain_variant_set_for_set(i) for i in SET_IDS]
    variant_sets_list = ga4gh_brca_page(variant_sets_list, int(page_size), int(page_token))
    if len(variant_sets_list) > page_size:
        variant_sets_list.pop()
        page_token = str(1 + int(page_token))
        response.next_page_token = page_token
    for sets in variant_sets_list:
        response.variant_sets.extend([sets])
    return JsonResponse(json_format.MessageToDict(response, True))

def obtain_variant_set_for_set(Set):
    variant_set = variants.VariantSet()
    variant_set.id = '{}-{}'.format(DATASET_ID, Set)
    variant_set.name = '{}-{}'.format(SETNAME, Set)
    variant_set.dataset_id = DATASET_ID
    variant_set.reference_set_id = '{}-{}'.format(REFERENCE_SET_BASE, Set)
    brca_meta(variant_set.metadata, DATASET_ID)
    return variant_set

def brca_meta(metadata, dataset_id):
    """Auxiliary function, generates metadata fields"""
    metadata_element = variants.VariantSetMetadata()
    for key in Variant._meta.get_fields():
        metadata_element.key = str(key.name)
        metadata_element.value = '-'
        metadata_element.id = '{}-{}'.format(dataset_id , str(key.name))
        metadata_element.type = key.get_internal_type()
        metadata_element.number = '-'
        metadata_element.description = "refer to ->{} in https://github.com/BD2KGenomics" \
                                       "/brca-website/blob/master/content/help_research.md".format(str(key.name))
        metadata.extend([metadata_element])
    return metadata


@require_http_methods(['GET'])
def get_variant_set(request, variant_set_id):
    """/variantsets/<set id> method"""
    if not variant_set_id:
        return HttpResponseBadRequest(
            json.dumps(ErrorMessages['variantSetId']),
            content_type='application/json')
    dataset, id_ = variant_set_id.split('-')

    if id_ in SET_IDS and dataset == 'brca':
        variant_set = variants.VariantSet()
        variant_set.id = '{}-{}'.format(dataset, id_)
        variant_set.name = '{}-{}'.format(SETNAME, id_)
        variant_set.dataset_id = DATASET_ID
        variant_set.reference_set_id = '{}-{}'.format(REFERENCE_SET_BASE, id_)
        brca_meta(variant_set.metadata, id_)
        resp = json_format.MessageToDict(variant_set, True)
        return JsonResponse(resp)
    else:
        return JsonResponse({'Invalid Set Id': variant_set_id}, status=404)


@require_http_methods(['POST'])
def search_datasets(request):
    """/datasets/search method request handler"""
    if not request.body:
        page_size = DEFAULT_PAGE_SIZE
        page_token = '0'
    else:
        try :
            request_dict = json_format.Parse(request.body, metadata_service.SearchDatasetsRequest())
        except Exception as e:
            return HttpResponseBadRequest(json.dumps({"message": e.message.replace("\"", "'")}), content_type='application/json')
        page_size = request_dict.page_size
        page_token = request_dict.page_token

    if not page_size or page_size == 0:
        page_size = DEFAULT_PAGE_SIZE
    if not page_token:
        page_token = '0'

    response = metadata_service.SearchDatasetsResponse()
    dataset = metadata.Dataset()
    dataset.name = SETNAME
    dataset.id = DATASET_ID
    dataset.description = 'Variants observed in brca-exchange project'
    # TODO uncomment when ga4gh client implements info field otherwise hardcoded values are placed
    #dataset.info[SETNAME].append("This set contains variants as stored and mantained by the brca-exchange project")
    list_of_sets = []
    list_of_sets.append(dataset)
    sets = ga4gh_brca_page(list_of_sets, int(page_size), int(page_token))
    response.next_page_token
    if len(sets) > page_size:
        page_token = str(1 + int(page_token))
        response.next_page_token = page_token
    ##############
    # TODO Block gets fixed when ga4gh client implements
    # info field otherwise hardcoded values are placed
    if sets:
        response.datasets.extend(sets)
    else:
        response.next_page_token = ' '
        response.datasets.extend([metadata.Dataset()])
    ##############
    return JsonResponse(json_format.MessageToDict(response, False))

@require_http_methods(['GET'])
def get_dataset(request, dataset_id):
    """/datasets/<dataset id> get dataset via id method"""
    if not dataset_id:
        return HttpResponseBadRequest(
            json.dumps(ErrorMessages['datasetId']),
            content_type='application/json')
    response = metadata_service.GetDatasetRequest()
    dataset = metadata.Dataset()
    response.dataset_id = dataset_id
    dataset.id = DATASET_ID
    dataset.name = SETNAME
    dataset.description = 'Variants observed in brca-exchange project'
    # Needs field for info, still not available from ga4gh client
    return JsonResponse(json_format.MessageToDict(dataset, False))

@require_http_methods(['GET', 'POST'])
def empty_variantset_id_catcher(request):
    """Error URL catcher methods"""
    return HttpResponseBadRequest(json.dumps(ErrorMessages['methodNotAllowed']),
                                  content_type='application/json',
                                  status=405)

@require_http_methods(['GET', 'POST'])
def empty_variant_id_catcher(request):
    return HttpResponseBadRequest(json.dumps(ErrorMessages['methodNotAllowed']),
                                  content_type='application/json',
                                  status=405)
@require_http_methods(['GET', 'POST'])
def empty_dataset_catcher(request):
    return HttpResponseBadRequest(json.dumps(ErrorMessages['methodNotAllowed']),
                                  content_type='application/json',
                                  status=405)

ErrorMessages = {'emptyBody' :{'status_code': 400, 'message' : 'Invalid request: empty request'},
                 'variantSetId' : {'status_code': 400, 'message': 'Invalid request: please provide a variantSetId'},
                 'referenceName': {'status_code': 400, 'message': 'Invalid request: please provide a referenceName, 13 or 17'},
                 'start': {'status_code' : 400, 'message': 'Invalid request: please provide a start position'},
                 'end' : {'status_code' :400, 'message': 'Invalid request: please provide an end position'},
                 'datasetId': {'status_code' : 400, 'message': 'Invalid request: please provide a datasetId'},
                 'variantId': {'status_code' : 400, 'message': 'Invalid request: please provide a variantId'},
                 'invalidPositions': {'status_code': 400, 'message': 'Invalid request: assure starting position is less than end'},
                 'notFoundId': {'status_code' : 404, 'message': 'Not found: the provided id is not supported'},
                 'methodNotAllowed': {'status_code': 405, 'message': 'Method is not supported: empty body request'}}

# The display name for the variant set.
SETNAME = 'brca-exchange-variants'

# The identifier for the dataset
DATASET_ID = 'brca'

# The string identify the reference set. Currently for display only.
REFERENCE_SET_BASE = 'Genomic-Coordinate'

# The name of the dataset for display
DATASET_NAME = 'brca-exchange'

# The list of reference genomes used to switch between variant sets.
SET_IDS = ['hg36', 'hg37', 'hg38']

# The list of reference names to be served
REFERENCE_NAMES = ['chr13', 'chr17', '13', '17']

# When no pagesize is specified pages of this length will be returned.
DEFAULT_PAGE_SIZE = 3

# Need to implement function that filters elements by increasing and decreasing integer values
