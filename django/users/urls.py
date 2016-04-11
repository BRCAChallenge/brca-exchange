from django.conf import settings
from django.conf.urls import url
from rest_framework_jwt.views import obtain_jwt_token

from . import views


urlpatterns = [
    url(r'^register/$', views.register, name="register"),
    url(r'^users/$', views.users, name="users"),
    url(r'^token-auth/$', obtain_jwt_token, name="token_auth"),
    url(r'^get/$', views.retrieve, name="retrieve"),
    url(r'^update/$', views.update, name="update"),
]
