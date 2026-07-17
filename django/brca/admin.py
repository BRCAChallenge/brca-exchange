from django.contrib import admin
from django.contrib.admin.sites import NotRegistered
from django.contrib.auth.models import Group
from django.contrib.sites.models import Site
from rest_framework.authtoken.models import Token

try:
    admin.site.unregister(Group)
except NotRegistered:
    pass
try:
    admin.site.unregister(Site)
except NotRegistered:
    pass
try:
    admin.site.unregister(Token)
except NotRegistered:
    pass

