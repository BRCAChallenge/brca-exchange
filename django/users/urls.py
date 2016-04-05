from django.conf import settings
from django.conf.urls import url, patterns

from . import views

url(r'^media/(?P<path>.*)$', 'django.views.static.serve', {'document_root': settings.MEDIA_ROOT}),

urlpatterns = [
    url(r'^register/$', views.register, name="register"),
    url(r'^login/$', views.login, name="login"),
    url(r'^logout/$', views.user_logout, name="logout"),
    url(r'^users/$', views.users, name="users"),

]
