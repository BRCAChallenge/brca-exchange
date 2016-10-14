from django.contrib import admin
from django.contrib.auth.models import Group
from django.contrib.sites.models import Site
from rest_framework.authtoken.models import Token

admin.site.unregister(Group)
admin.site.unregister(Site)
admin.site.unregister(Token)

