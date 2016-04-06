from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^register/$', views.register, name="register"),
    url(r'^users/$', views.users, name="users"),
    url(r'^token-auth/$', views.token_auth, name="token_auth"),
]
