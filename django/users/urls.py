from django.conf.urls import url
from rest_framework_jwt.views import obtain_jwt_token

from . import views

urlpatterns = [
    url(r'^register/$', views.register, name="register"),
    url(r'^mailinglist/$', views.mailinglist, name="mailinglist"),
    url(r'^users/$', views.users, name="users"),
    url(r'^user_locations/$', views.user_locations, name="user_locations"),
    url(r'^token-auth/$', obtain_jwt_token, name="token_auth"),
    url(r'^get/$', views.retrieve, name="retrieve"),
    url(r'^update/$', views.update, name="update"),
    url(r'^password_reset/$', views.password_reset, name="password_reset"),
    url(r'^confirm/(?P<activation_key>\w+)/$', views.confirm, name="confirm"),
    url(r'^resend-activation/$', views.resend_activation, name="resend_activation"),
    url(r'^update_password/(?P<password_reset_token>\w+)/$', views.update_password, name="update_password"),
    url(r'^check_password_token/(?P<password_reset_token>\w+)/$', views.check_password_token, name="check_password_token"),
]
