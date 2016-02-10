import json

from django.core import serializers
from django.http import JsonResponse

from .models import Variant


def index(request):
    # data = json.load(open('/home/pete/work/BRCA/brca-website/brca.json','r'))
    querySet = Variant.objects.all()
    json_data = serializers.serialize('json', querySet)

    data = json.loads(json_data)

    header = data[0]['fields'].keys()
    rows = [row['fields'].values() for row in data]

    response_data = {'header': header, 'rows': rows}

    response = JsonResponse(response_data)

    response['Access-Control-Allow-Origin'] = '*'

    return response
