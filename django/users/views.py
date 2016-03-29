from django.contrib import auth
from django.contrib.auth import logout
from django.db import IntegrityError
from django.http import JsonResponse

from .models import MyUser


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

            return JsonResponse({'success': True})

    else:
        return JsonResponse({'success': False})


def user_logout(request):
    logout(request)
    return JsonResponse({'success': True})


def register(request):
    email = request.POST.get('email', '')
    password = request.POST.get('password', '')

    title = request.POST.get('title', '')
    affiliation = request.POST.get('affiliation', '')
    institution = request.POST.get('institution', '')
    city = request.POST.get('city', '')
    state = request.POST.get('state', '')
    country = request.POST.get('country', '')
    phone_number = request.POST.get('phone_number', '')
    hide_number = request.POST.get('hide_number', True)
    hide_email = request.POST.get('hide_email', True)

    try:
        MyUser.objects.create_user(email, password, title, affiliation, institution, city, state, country, phone_number,
                                   hide_number, hide_email)
    except IntegrityError:
        return JsonResponse({'success': False})

    return JsonResponse({'success': True})
