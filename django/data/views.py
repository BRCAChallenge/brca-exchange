import re
from operator import __or__
import tempfile
import os

from django.db import connection
from django.db.models import Q
from django.http import JsonResponse, HttpResponse

from .models import Variant


def sanitise_term(term):
    # Escape all non alphanumeric characters
    term = re.escape(term)
    # Enable prefix search
    term += ":*"
    return term


def index(request):
    order_by = request.GET.get('order_by')
    direction = request.GET.get('direction')
    page_size = int(request.GET.get('page_size', '0'))
    page_num = int(request.GET.get('page_num', '0'))
    search_term = request.GET.get('search_term')
    format = request.GET.get('format')
    source = request.GET.getlist('source')
    filters = request.GET.getlist('filter')
    filter_values = request.GET.getlist('filterValue')

    query, count = build_query(direction, filter_values, filters, order_by, search_term, source, page_size, page_num)
    # The query for producing the CSV needs to have string literals inside single quotes
    query_csv, _ = build_query(direction, filter_values, filters, order_by, search_term, source, page_size, page_num,
                               quotes="\'")

    if format == 'csv':
        cursor = connection.cursor()
        with tempfile.NamedTemporaryFile() as f:
            os.chmod(f.name, 0606)
            cursor.execute("COPY ({}) TO '{}' WITH DELIMITER ',' CSV HEADER".format(query_csv.query, f.name))

            response = HttpResponse(f.read(), content_type='text/csv')
            response['Content-Disposition'] = 'attachment;filename="variants.csv"'
            return response

    elif format == 'json':
        # call list() now to evaluate the query
        response = JsonResponse({'count': count, 'data': list(query.values())})
        response['Access-Control-Allow-Origin'] = '*'
        return response


def build_query(direction, filterValues, filters, order_by, search_term, source, page_size, page_num,
                quotes=""):
    query = Variant.objects

    # if there are multiple filters the row must match all the filters
    if filters:
        for column, value in zip(filters, filterValues):
            if column == 'id':
                query = query.filter(**{column:value})
            else:
                query = query.extra(
                    where=["\"{0}\" LIKE %s".format(column)],
                    params=["{0}{1}%{0}".format(quotes, value)]
                )

    # search using the tsvector column which represents our document made of all the columns
    if search_term:
        query = query.extra(
            where=["variant.fts_document @@ to_tsquery('simple', %s)"],
            params=["{0}{1}{0}".format(quotes, sanitise_term(search_term))]
        )

    # if there are multiple sources given then OR them:
    # the row must match in at least one column
    if source:
        query_list = (Q(**{column: True}) for column in source)
        query = query.filter(reduce(__or__, query_list))
    if order_by:
        # special case for HGVS columns
        if order_by in ('HGVS_cDNA','HGVS_protein'):
            order_by = 'Genomic_Coordinate'
        if direction == 'descending':
            order_by = '-' + order_by
        query = query.order_by(order_by)

    # count the number of rows now before paginating
    count = query.count()

    if page_size:
        start = page_size * page_num
        end = start + page_size
        query = query[start:end]
    return query, count


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
