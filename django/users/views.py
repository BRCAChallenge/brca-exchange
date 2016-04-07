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
        user = request.user
        query = MyUser.objects.filter(email=user)
        data = list(query.values())[0]
        data["password"] = ""
        response = JsonResponse({'user': data})
        return response


class UpdateUser(APIView):
    permission_classes = (IsAuthenticated,)
    authentication_classes = (JSONWebTokenAuthentication,)

    def update(self, request):
        user = request.user
        user = MyUser.objects.filter(email=user)
        # user.update(...)

        response = JsonResponse({'success': True})
        return response


def token_auth(request):
    response = obtain_jwt_token(request)
    return response


def register(request):
    fields = user_fields(request)
    image = None
    if request.FILES:
        image = request.FILES["image"]
    # Check the CAPTCHA
    try:
        captcha = request.POST.get('captcha')
        post_data = {'secret': settings.CAPTCHA_SECRET,
                     'response': captcha}
        response = requests.post('https://www.google.com/recaptcha/api/siteverify', data=post_data)
        content = json.loads(response.content)
        response = {'success': content['success']}
    except HTTPError:
        response = {'success': False}

    # Create the user
    try:
        created_user = MyUser.objects.create_user(**fields)
        # Save the image under the user's id
        if image is not None:
            save_picture(created_user.id, image)

    except IntegrityError:
        response = {'success': False}

    response = JsonResponse(response)
    return response


def user_fields(request):
    image = None
    if request.FILES:
        image = request.FILES["image"]
    has_image = (image is not None)
    email = request.POST.get('email', '')
    password = request.POST.get('password', '')
    firstName = request.POST.get('firstName', '')
    lastName = request.POST.get('lastName', '')
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

    return {'affiliation': affiliation, 'city': city, 'comment': comment, 'country': country,
            'email': email, 'firstName': firstName, 'has_image': has_image, 'hide_email': hide_email,
            'hide_number': hide_number, 'include_me': include_me, 'institution': institution,
            'lastName': lastName, 'password': password, 'phone_number': phone_number, 'state': state, 'title': title}


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
