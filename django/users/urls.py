from django.urls import re_path
from rest_framework_simplejwt.views import TokenObtainPairView

from . import views

urlpatterns = [
    re_path(r'^register/$', views.register, name="register"),
    re_path(r'^users/$', views.users, name="users"),
    re_path(r'^user_locations/$', views.user_locations, name="user_locations"),
    re_path(r'^token-auth/$', TokenObtainPairView.as_view(), name="token_auth"),
    re_path(r'^get/$', views.retrieve, name="retrieve"),
    re_path(r'^update/$', views.update, name="update"),
    re_path(r'^password_reset/$', views.password_reset, name="password_reset"),
    re_path(r'^confirm/(?P<activation_key>\w+)/$', views.confirm, name="confirm"),
    re_path(r'^resend-activation/$', views.resend_activation, name="resend_activation"),
    re_path(r'^update_password/(?P<password_reset_token>\w+)/$', views.update_password, name="update_password"),
    re_path(r'^check_password_token/(?P<password_reset_token>\w+)/$', views.check_password_token, name="check_password_token"),
]
