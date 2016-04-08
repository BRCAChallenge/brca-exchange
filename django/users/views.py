import datetime
import hashlib
import json
import os
import random
from urllib2 import HTTPError

import requests
from django.core.mail import send_mail
from django.db import IntegrityError
from django.http import HttpResponse
from django.http import JsonResponse
from django.shortcuts import get_object_or_404
from django.views.decorators.cache import never_cache
from rest_framework.decorators import api_view, permission_classes, authentication_classes
from rest_framework.permissions import IsAuthenticated
from rest_framework_jwt.authentication import JSONWebTokenAuthentication

from brca import settings
from .models import MyUser


@api_view(['GET'])
@permission_classes((IsAuthenticated,))
@authentication_classes((JSONWebTokenAuthentication,))
def retrieve(request):
    user = request.user
    query = MyUser.objects.filter(email=user)
    data = list(query.values())[0]
    data["password"] = ''
    response = JsonResponse({'user': data})
    return response


@never_cache
@api_view(['POST'])
@permission_classes((IsAuthenticated,))
@authentication_classes((JSONWebTokenAuthentication,))
def update(request):
    print(request)
    user = MyUser.objects.filter(email=request.user)

    fields = user_fields(request)
    delete_image = request.POST.get('deleteImage', "false") == "true"
    fields['has_image'] = (user[0].has_image or fields['has_image']) and not delete_image

    del fields['email']
    if fields['password'] == '':
        del fields['password']

    try:
        user.update(**fields)
        if request.FILES:
            image = request.FILES["image"]
            if image is not None:
                save_picture(user[0].id, image)
    except Exception, e:
        return JsonResponse({'success': False, 'error': str(e)})

    return JsonResponse({'success': True})


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

        created_user.is_active = False
        created_user.save()

        # Send email with activation key
        salt = hashlib.sha1(str(random.random())).hexdigest()[:5]
        activation_key = hashlib.sha1(salt + created_user.email).hexdigest()
        key_expires = datetime.datetime.today() + datetime.timedelta(2)

        created_user.activation_key = activation_key
        created_user.key_expires = key_expires
        created_user.save()

        email_subject = 'Account confirmation'
        email_body = "Hey %s, thanks for signing up. To activate your account, click this link within 48 hours http://127.0.0.1:8000/accounts/confirm/%s" % (
            created_user.firstName, activation_key)
        send_mail(email_subject, email_body, 'noreply@brcaexchange.org', [created_user.email], fail_silently=False)

    except IntegrityError:
        response = {'success': False}

    return JsonResponse(response)


def confirm(request, activation_key):
    user = get_object_or_404(MyUser, activation_key=activation_key)
    user.is_active = True
    user.save()
    return HttpResponse("Thanks for confirming your email")


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
    phone_number = request.POST.get('phone_number', '')
    comment = request.POST.get('comment', '')
    include_me = (request.POST.get('include_me', "true") == "true")
    hide_number = (request.POST.get('hide_number', "true") == "true")
    hide_email = (request.POST.get('hide_email', "true") == "true")

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

    count = query.count()

    start = page_num * page_size
    end = start + page_size
    data = list(query[start:end].values())

    for user in data:
        if user['hide_email']:
            user['email'] = ""
        if user['hide_number']:
            user['phone_number'] = ""

    response = JsonResponse({'data': data, 'count': count})
    return response
