from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^$', views.index, name="index"),
    url(r'^releases', views.releases, name="releases"),
    url(r'^variant/$', views.variant, name="variant"),
    url(r'^suggestions/$', views.autocomplete)
]
