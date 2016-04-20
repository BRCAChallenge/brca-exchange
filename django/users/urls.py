from django.conf.urls import url
from rest_framework_jwt.views import obtain_jwt_token

from . import views

urlpatterns = [
    url(r'^register/$', views.register, name="register"),
    url(r'^users/$', views.users, name="users"),
    url(r'^token-auth/$', obtain_jwt_token, name="token_auth"),
    url(r'^get/$', views.retrieve, name="retrieve"),
    url(r'^update/$', views.update, name="update"),
    url(r'^password_reset/$', views.password_reset, name="password_reset"),
    url(r'^confirm/(?P<activation_key>\w+)/$', views.confirm, name="confirm"),
    url(r'^update_password/(?P<password_reset_token>\w+)/$', views.update_password, name="update_password"),
]
