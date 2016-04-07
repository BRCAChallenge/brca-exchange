import json
import os
from urllib2 import HTTPError

import requests
from django.db import IntegrityError
from django.http import JsonResponse
from rest_framework.permissions import IsAuthenticated
from rest_framework.response import Response
from rest_framework.views import APIView
from rest_framework_jwt.authentication import JSONWebTokenAuthentication
from rest_framework_jwt.views import obtain_jwt_token

from brca import settings
from .models import MyUser


class RetrieveUser(APIView):
    permission_classes = (IsAuthenticated,)
    authentication_classes = (JSONWebTokenAuthentication,)

    def get(self, request):
        data = {
            'token': str(request.auth)
        }
        response = Response(data)
        return response


def token_auth(request):
    response = obtain_jwt_token(request)
    return response


def register(request):
    image = None
    if request.FILES:
        image = request.FILES["image"]
    has_image = (image is not None)

    email = request.POST.get('email', '')
    password = request.POST.get('password', '')

    first_name = request.POST.get('firstName', '')
    last_name = request.POST.get('lastName', '')
    title = request.POST.get('title', '')
    affiliation = request.POST.get('affiliation', '')
    institution = request.POST.get('institution', '')
    city = request.POST.get('city', '')
    state = request.POST.get('state', '')
    country = request.POST.get('country', '')
    phone_number = request.POST.get('phoneNumber', '')
    comment = request.POST.get('comment', '')
    include_me = (request.POST.get('includeMe', "true") == "true")
    hide_number = (request.POST.get('hideNumber', "true") == "true")
    hide_email = (request.POST.get('hideEmail', "true") == "true")
    captcha = request.POST.get('captcha')

    response = {'success': True}

    # Check the CAPTCHA
    try:
        post_data = {'secret': settings.CAPTCHA_SECRET,
                     'response': captcha}
        response = requests.post('https://www.google.com/recaptcha/api/siteverify', data=post_data)
        content = json.loads(response.content)
        response = {'success': content['success']}
    except HTTPError:
        response = {'success': False}

    # Create the user
    try:
        created_user = MyUser.objects.create_user(email=email,
                                                  password=password,
                                                  firstName=first_name,
                                                  lastName=last_name,
                                                  title=title,
                                                  affiliation=affiliation,
                                                  institution=institution,
                                                  city=city,
                                                  state=state,
                                                  comment=comment,
                                                  country=country,
                                                  phone_number=phone_number,
                                                  include_me=include_me,
                                                  hide_number=hide_number,
                                                  hide_email=hide_email,
                                                  has_image=has_image)
        # Save the image under the user's id
        if image is not None:
            save_picture(created_user.id, image)

    except IntegrityError:
        response = {'success': False}

    response = JsonResponse(response)
    return response


def save_picture(filename, image):
    path = os.path.join(settings.MEDIA_ROOT, str(filename))
    fd = open(path, 'wb')
    for chunk in image.chunks():
        fd.write(chunk)
    fd.close()


def users(request):
    page_num = int(request.GET.get('page_num', '0'))
    page_size = int(request.GET.get('page_size', '0'))

    query = MyUser.objects.filter(include_me=True)

    start = page_num * page_size
    end = start + page_size
    data = list(query[start:end].values())

    for user in data:
        if user['hide_email']:
            user['email'] = ""
        if user['hide_number']:
            user['phone_number'] = ""

    response = JsonResponse({'data': data})
    return response
