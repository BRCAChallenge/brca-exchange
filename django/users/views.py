import requests, json

from django.contrib import auth
from django.contrib.auth import logout
from django.db import IntegrityError
from django.http import JsonResponse

from .models import MyUser
from .config import captcha_secret


def login(request):
    # here you get the post request username and password
    username = request.POST.get('username', '')
    password = request.POST.get('password', '')

    # authentication of the user, to check if it's active or None
    user = auth.authenticate(username=username, password=password)

    if user is not None:
        if user.is_active:
            # this is where the user login actually happens, before this the user
            # is not logged in.
            auth.login(request, user)

            response = JsonResponse({'success': True})
            response["Access-Control-Allow-Origin"] = "*"
            return response

    else:
        response = JsonResponse({'success': False})
        response["Access-Control-Allow-Origin"] = "*"
        return response


def user_logout(request):
    logout(request)
    response = JsonResponse({'success': True})
    response["Access-Control-Allow-Origin"] = "*"
    return response


def register(request):
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
    include_me = request.POST.get('includeMe', True)
    hide_number = request.POST.get('hideNumber', True)
    hide_email = request.POST.get('hideEmail', True)
    captcha = request.POST.get('captcha')

    response = {'success': True}

    # Check the CAPTCHA
    try:
        post_data = {'secret': captcha_secret,
                     'response': captcha}
        response = requests.post('https://www.google.com/recaptcha/api/siteverify', data=post_data)
        content = json.loads(response.content)
        response = {'success': content['success']}
    except HTTPError:
        response = {'success': False}

    try:
        MyUser.objects.create_user(email, password, first_name, last_name, title, affiliation, institution, city, state,
                                   comment, country, phone_number, include_me, hide_number, hide_email)
    except IntegrityError:
        response = {'success': False}

    response = JsonResponse(response)
    response["Access-Control-Allow-Origin"] = "*"
    return response

def users(request):
    page_num = int(request.GET.get('page_num','0'))
    page_size = int(request.GET.get('page_size','0'))

    start = page_num * page_size
    end = start + page_size

    page = MyUser.objects.all()[start:end]

    response = JsonResponse({'data':list(page.values())})
    response["Access-Control-Allow-Origin"] = "*"
    
    return response
