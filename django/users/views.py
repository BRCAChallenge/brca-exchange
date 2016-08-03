import datetime
import hashlib
import json
import os
import random
from urllib2 import HTTPError

import requests
from django.core.mail import EmailMultiAlternatives
from django.db import IntegrityError
from django.db.models import Q
from django.http import JsonResponse
from django.template import Context
from django.template.loader import get_template
from django.utils import timezone
from django.views.decorators.cache import never_cache
from rest_framework.decorators import api_view, permission_classes, authentication_classes
from rest_framework.permissions import IsAuthenticated
from rest_framework_jwt.authentication import JSONWebTokenAuthentication

from brca import settings, site_settings
from .models import MyUser, MailingListEmail


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
        return JsonResponse({'success': False, 'error': 'Wrong CAPTCHA'})

    # Create the user
    try:
        created_user = MyUser.objects.create_user(**fields)
        # Save the image under the user's id
        if image is not None:
            save_picture(created_user.id, image)

        created_user.is_active = False
        created_user.save()

        # Create and save activation key
        salt = hashlib.sha1(str(random.random())).hexdigest()[:5]
        activation_key = hashlib.sha1(salt + created_user.email).hexdigest()

        created_user.activation_key = activation_key
        created_user.save()

        # Send activation email
        url = "{0}confirm/{1}".format(site_settings.URL_FRONTEND, activation_key)
        plaintext_email = get_template(os.path.join(settings.BASE_DIR, 'users', 'templates', 'registration_email.txt'))
        html_email = get_template(os.path.join(settings.BASE_DIR, 'users', 'templates', 'registration_email.html'))

        d = Context({'firstname': created_user.firstName, 'url': url})

        subject, from_email, to = 'BRCAExchange account confirmation', 'noreply@brcaexchange.org', created_user.email
        text_content = plaintext_email.render(d)
        html_content = html_email.render(d)
        msg = EmailMultiAlternatives(subject, text_content, from_email, [to])
        msg.attach_alternative(html_content, "text/html")
        msg.send()

    except IntegrityError:
        return JsonResponse({'success': False, 'error': 'This email address already exists'})

    return JsonResponse({'success': True})

def resend_activation(request):
    email = request.POST.get('email', '')
    user = MyUser.objects.filter(email = email)
    if not user:
        response = JsonResponse({'success': False, 'error': 'Email not found'})
        response['Access-Control-Allow-Origin'] = '*'
        return response;
    user = user[0]

    # resend activation email
    url = "{0}confirm/{1}".format(site_settings.URL_FRONTEND, user.activation_key)
    plaintext_email = get_template(os.path.join(settings.BASE_DIR, 'users', 'templates', 'registration_email.txt'))
    html_email = get_template(os.path.join(settings.BASE_DIR, 'users', 'templates', 'registration_email.html'))

    d = Context({'firstname': user.firstName, 'url': url})

    subject, from_email, to = 'BRCAExchange account confirmation', 'noreply@brcaexchange.org', user.email
    text_content = plaintext_email.render(d)
    html_content = html_email.render(d)
    msg = EmailMultiAlternatives(subject, text_content, from_email, [to])
    msg.attach_alternative(html_content, "text/html")
    msg.send()

    return JsonResponse({'success': True})

def confirm(request, activation_key):
    user = MyUser.objects.filter(activation_key=activation_key)
    if not user:
        response = JsonResponse({'success': False, 'error': 'Invalid activation key'})
        response['Access-Control-Allow-Origin'] = '*'
        return response
    user = user[0]
    user.is_active = True
    user.save()
    response = JsonResponse({'success': True})
    response['Access-Control-Allow-Origin'] = '*'
    return response


def password_reset(request):
    email = request.POST.get('email', '')
    user = MyUser.objects.filter(email=email)
    if not user:
        response = JsonResponse({'success': False, 'error': 'Email not found'})
        response['Access-Control-Allow-Origin'] = '*'
        return response
    user = user[0]

    # Create and save password reset token
    salt = hashlib.sha1(str(random.random())).hexdigest()[:5]
    password_reset_token = hashlib.sha1(salt + user.email).hexdigest()

    user.password_reset_token = password_reset_token
    email_duration_days = settings.PASSWORD_RESET_LINK_DURATION
    password_token_expires = timezone.now() + datetime.timedelta(email_duration_days)
    user.password_token_expires = password_token_expires
    user.save()

    # Send password reset email
    url = "{0}reset/{1}".format(site_settings.URL_FRONTEND, password_reset_token)
    plaintext_email = get_template(os.path.join(settings.BASE_DIR, 'users', 'templates', 'password_reset_email.txt'))
    html_email = get_template(os.path.join(settings.BASE_DIR, 'users', 'templates', 'password_reset_email.html'))
    d = Context({'firstname': user.firstName, 'url': url, 'hours': email_duration_days * 24})

    subject, from_email, to = 'BRCAExchange password reset', 'noreply@brcaexchange.org', user.email
    text_content = plaintext_email.render(d)
    html_content = html_email.render(d)
    msg = EmailMultiAlternatives(subject, text_content, from_email, [to])
    msg.attach_alternative(html_content, "text/html")
    msg.send()

    response = JsonResponse({'success': True})
    response['Access-Control-Allow-Origin'] = '*'
    return response


def check_password_token(request, password_reset_token):
    user = MyUser.objects.filter(password_reset_token=password_reset_token)
    if not user or user[0].password_token_expires < timezone.now():
        response = JsonResponse({'invalid_token': True})
        response['Access-Control-Allow-Origin'] = '*'
        return response

    response = JsonResponse({'invalid_token': False})
    response['Access-Control-Allow-Origin'] = '*'
    return response


def update_password(request, password_reset_token):
    user = MyUser.objects.filter(password_reset_token=password_reset_token)
    if not user:
        response = JsonResponse({'success': False, 'invalid_token': True})
        response['Access-Control-Allow-Origin'] = '*'
        return response
    user = user[0]

    if user.password_token_expires < timezone.now():
        response = JsonResponse({'success': False, 'invalid_token': True})
        response['Access-Control-Allow-Origin'] = '*'
        return response

    password = request.POST.get('password', '')
    user.set_password(password)
    user.password_token_expires = timezone.now()
    user.save()
    response = JsonResponse({'success': True})
    response['Access-Control-Allow-Origin'] = '*'
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
    phone_number = request.POST.get('phone_number', '')
    comment = request.POST.get('comment', '')
    include_me = (request.POST.get('include_me', "true") == "true")
    email_me = (request.POST.get('email_me', "true") == "true")
    hide_number = (request.POST.get('hide_number', "true") == "true")
    hide_email = (request.POST.get('hide_email', "true") == "true")

    return {'affiliation': affiliation, 'city': city, 'comment': comment, 'country': country,
            'email': email, 'firstName': firstName, 'has_image': has_image, 'hide_email': hide_email,
            'hide_number': hide_number, 'include_me': include_me, 'email_me': email_me, 'institution': institution,
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
    search = request.GET.get('search', '')

    query = MyUser.objects.filter(include_me=True).filter(is_approved=True)
    search_query = Q()
    for term in search.split():
        search_query &= Q(firstName__icontains=term) | Q(lastName__icontains=term) | Q(institution__icontains=term) | Q(city__icontains=term) | Q(state__icontains=term) | Q(country__icontains=term)
    query = query.filter(search_query)

    whitelist = ['id','email','firstName','lastName','title','affiliation','institution','city','state','country','phone_number','hide_number','hide_email','include_me','is_active','is_admin','comment','has_image','is_approved','email_me']

    count = query.count()

    start = page_num * page_size
    end = start + page_size
    data = list(query[start:end].values(*whitelist))


    for user in data:
        if user['hide_email']:
            user['email'] = ""
        if user['hide_number']:
            user['phone_number'] = ""

    response = JsonResponse({'data': data, 'count': count})
    return response

def user_locations(request):
    query = MyUser.objects.filter(include_me=True).filter(is_approved=True)
    fields = ['id', 'firstName', 'lastName', 'title', 'institution', 'city', 'state', 'country', 'has_image']
    response = JsonResponse({'data': list(query.values(*fields))})
    return response

def mailinglist(request):
    email = request.POST.get('email')

    # Check the CAPTCHA
    try:
        captcha = request.POST.get('captcha')
        post_data = {'secret': settings.CAPTCHA_SECRET,
                     'response': captcha}
        response = requests.post('https://www.google.com/recaptcha/api/siteverify', data=post_data)
        content = json.loads(response.content)
        response = {'success': content['success']}
    except HTTPError:
        return JsonResponse({'success': False, 'error': 'Wrong CAPTCHA'})

    try:
        entry = MailingListEmail.objects.create(email = request.POST.get('email'))
        entry.save()
    except IntegrityError:
        return JsonResponse({'success': False, 'error': 'This email address already exists'})

    # Send confirmation email
    plaintext_email = get_template(os.path.join(settings.BASE_DIR, 'users', 'templates', 'mailinglist_email.txt'))
    html_email = get_template(os.path.join(settings.BASE_DIR, 'users', 'templates', 'mailinglist_email.html'))

    subject, from_email, to = 'BRCA Exchange Mailing List', 'noreply@brcaexchange.org', email
    text_content = plaintext_email.render()
    html_content = html_email.render()
    msg = EmailMultiAlternatives(subject, text_content, from_email, [to])
    msg.attach_alternative(html_content, "text/html")
    msg.send()

    return JsonResponse({'success': True}) 
