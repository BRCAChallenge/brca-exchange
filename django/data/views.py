import json

from operator import __or__

from django.core import serializers
from django.http import JsonResponse
from django.db.models import Q

from .models import Variant

def index(request):
    
    order_by =       request.GET.get('order_by')
    direction =      request.GET.get('direction')
    page_size =  int(request.GET.get('page_size'))
    page_num =   int(request.GET.get('page_num'))
    search_term =    request.GET.get('search_term')
    search_columns = request.GET.getlist('search_column')
    source =         request.GET.getlist('source')
    filters =        request.GET.getlist('filter')
    filterValues =   request.GET.getlist('filterValue')

    query = Variant.objects.values()

    # if there are multiple filters given then AND them:
    # the row must match all the filters
    if filters:
        query = query.filter(**dict(zip(filters,filterValues)))

    # if there are multiple search columns given then OR them:
    # the row must match in at least one column
    if search_term:
        query_list = (Q(**{column+'__icontains':search_term}) for column in search_columns)
        query = query.filter(reduce(__or__,query_list))

    # if there are multiple sources given then OR them:
    # the row must match in at least one column
    if source:
        query_list = (Q(**{column:True}) for column in source)
        query = query.filter(reduce(__or__,query_list))

    # count the number of rows now before paginating
    count = query.count()

    if order_by:
        if direction == 'descending':
            order_by = '-' + order_by
        query = query.order_by(order_by)

    if page_size:
        start = page_size * page_num
        end = start + page_size
        query = query[start:end]

    # call list() now to evaluate the query
    response = JsonResponse({'count':count, 'data':list(query)})

    response['Access-Control-Allow-Origin'] = '*'

    return response
