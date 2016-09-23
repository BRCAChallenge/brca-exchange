from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^$', views.index, name="index"),
    url(r'^releases', views.releases, name="releases"),
    url(r'^suggestions/$', views.autocomplete)
]
