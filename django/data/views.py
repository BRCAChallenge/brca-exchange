import json

from django.core import serializers
from django.http import JsonResponse

from .models import Variant


def index(request):

    # Convert django's model representation to the format
    # expected by the frontend

    header = map(lambda field: field.name, Variant._meta.get_fields())
    rows = list(Variant.objects.values_list())

    data = {'header': header, 'rows': rows}

    response = JsonResponse(data)

    response['Access-Control-Allow-Origin'] = '*'

    return response
