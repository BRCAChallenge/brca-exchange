"""brca URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/5.1/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.conf import settings
from django.urls import path, include, re_path
from django.contrib import admin
from django.views.static import serve as static_serve
from rest_framework_simplejwt.views import (
    TokenObtainPairView,
    TokenRefreshView,
    )
urlpatterns = [
    path("api/token/", TokenObtainPairView.as_view(), name="token_obtain_pair"),
    path("api/token/refresh/", TokenRefreshView.as_view(), name="token_refresh"),
    re_path(r'^data/', include('data.urls')),
    re_path(r'^accounts/', include('users.urls')),
    re_path(r'^admin/', admin.site.urls),

    re_path(r'^site_media/media/(?P<path>.*)$', static_serve,
        {'document_root': settings.MEDIA_ROOT, 'show_indexes': True}),
    re_path(r'^downloads/(?P<path>.*)$', static_serve,
        {'document_root': settings.DOWNLOADS_ROOT, 'show_indexes': True}),
    re_path(r'^static/(?P<path>.*)$', static_serve,
        {'document_root': settings.STATIC_ROOT, 'show_indexes': True})
]
