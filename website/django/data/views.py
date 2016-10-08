import os
import re
import tempfile
from operator import __or__

from django.db import connection
from django.db.models import Q
from django.forms.models import model_to_dict
from django.http import JsonResponse, HttpResponse
from django.views.decorators.gzip import gzip_page

from .models import Variant, DataRelease, ChangeType

def releases(request):
    releases = DataRelease.objects.values().all()
    change_types = {x['name']:x['id'] for x in ChangeType.objects.values().all()}
    for release in releases:
        variants = Variant.objects.filter(Data_Release_id = release['id'])
        release['variants_added'] = variants.filter(Change_Type_id = change_types['new']).count()
        release['variants_classified'] = variants.filter(
            Q(Change_Type_id = change_types['changed_classification']) | Q(Change_Type_id = change_types['added_classification'])).count()
        release['variants_modified'] = variants.filter(
            Q(Change_Type_id = change_types['added_information']) | Q(Change_Type_id = change_types['changed_information'])).count()
        release['variants_deleted'] = variants.filter(Change_Type_id = change_types['deleted']).count()
    response = JsonResponse(list(releases), safe=False)
    response['Access-Control-Allow-Origin'] = '*'
    return response

def variant(request):
    variant_id = int(request.GET.get('variant_id'))

    variant = Variant.objects.get(id = variant_id)
    key = variant.Genomic_Coordinate_hg38

    query = Variant.objects.filter(Genomic_Coordinate_hg38 = key).order_by('-Data_Release_id').select_related('Data_Release')

    variant_versions = map(variant_to_dict, query)
    response = JsonResponse({"data": variant_versions})
    response['Access-Control-Allow-Origin'] = '*'
    return response

def variant_to_dict(variant_object):
    variant_dict = model_to_dict(variant_object)
    variant_dict["Data_Release"] = model_to_dict(variant_object.Data_Release)
    variant_dict["Data_Release"]["date_released"] = variant_object.Data_Release.date_released
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

    if release:
        query = Variant.objects.filter(Data_Release_id = int(release))
    else:
        latest = Variant.objects.distinct('Genomic_Coordinate_hg38').order_by('Genomic_Coordinate_hg38', '-Data_Release_id')
        query = Variant.objects.filter(id__in = latest)

    if format == 'csv':
        quotes = '\''
    else:
        quotes = ''

    if include or exclude:
        query = apply_sources(query, include, exclude)

    if filters:
        query = apply_filters(query, filter_values, filters, quotes=quotes)

    if search_term:
        query = apply_search(query, search_term, quotes=quotes)

    if order_by:
        query = apply_order(query, order_by, direction)

    if format == 'csv':

        cursor = connection.cursor()
        with tempfile.NamedTemporaryFile() as f:
            os.chmod(f.name, 0606)
            cursor.execute("COPY ({}) TO '{}' WITH DELIMITER ',' CSV HEADER".format(query.query, f.name))

            response = HttpResponse(f.read(), content_type='text/csv')
            response['Content-Disposition'] = 'attachment;filename="variants.csv"'
            return response

    elif format == 'json':

        count = query.count()

        if search_term:
            # Number of synonym matches = total matches minus matches on "normal" columns
            synonyms = count - apply_search(query, search_term, search_column='fts_standard').count()
        else:
            synonyms = 0

        query = select_page(query, page_size, page_num)

        # call list() now to evaluate the query
        response = JsonResponse({'count': count, 'synonyms': synonyms, 'data': list(query.values(*column))})
        response['Access-Control-Allow-Origin'] = '*'
        return response


def apply_sources(query, include, exclude):
    # if there are multiple sources given then OR them:
    # the row must match in at least one column
    include_list = (Q(**{column: True}) for column in include)
    exclude_dict = {exclusion: False for exclusion in exclude}
    
    return query.filter(reduce(__or__, include_list)).filter(**exclude_dict)


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


def apply_search(query, search_term, search_column='fts_document', quotes=''):
    # search using the tsvector column which represents our document made of all the columns
    where_clause = "variant.{} @@ to_tsquery('simple', %s)".format(search_column)
    parameter = quotes + sanitise_term(search_term) + quotes
    return query.extra(
        where=[where_clause],
        params=[parameter]
    )


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
    term = request.GET.get('term')
    limit = int(request.GET.get('limit', 10))

    cursor = connection.cursor()

    cursor.execute(
        """SELECT word FROM words
        WHERE word LIKE %s
        AND char_length(word) >= 3
        ORDER BY word""",
        ["%s%%" % term])

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
