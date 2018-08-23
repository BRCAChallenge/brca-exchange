"""brca URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.9/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.conf.urls import url, include
    2. Add a URL to urlpatterns:  url(r'^blog/', include('blog.urls'))
"""
from django.conf import settings
from django.conf.urls import include
from django.conf.urls import url, patterns
from django.contrib import admin
from django.views.static import serve as static_serve

urlpatterns = [
    url(r'^data/', include('data.urls')),
    url(r'^accounts/', include('users.urls')),
    url(r'^admin/', admin.site.urls),

    url(r'^site_media/media/(?P<path>.*)$', static_serve,
        {'document_root': settings.MEDIA_ROOT, 'show_indexes': True}),
    url(r'^downloads/(?P<path>.*)$', static_serve,
        {'document_root': settings.DOWNLOADS_ROOT, 'show_indexes': True}),
    url(r'^static/(?P<path>.*)$', static_serve,
        {'document_root': settings.STATIC_ROOT, 'show_indexes': True})
]
