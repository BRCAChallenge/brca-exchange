from django.contrib import admin
from django.db import models
from django.forms import TextInput
from .models import MyUser
from brca.settings import MEDIA_URL

def image(obj):
    if(obj.has_image):
        return '<img width="200" src="%s%s"/>' % (MEDIA_URL, obj.id)
    else:
        return ''
image.allow_tags = True

class MyUserAdmin(admin.ModelAdmin):

    list_display = ('email', 'is_admin', 'firstName', 'lastName', 'comment', 'is_approved', image)
    exclude = ('password',)
    ordering = ('is_approved',)
    search_fields = ('firstName', 'lastName', 'email')
    formfield_overrides = {
        models.TextField: {'widget': TextInput(attrs={'size':'20'})}
    }

    actions = ('approve',)

    def approve(self, request, queryset):
        queryset.update(is_approved=True)

admin.site.register(MyUser, MyUserAdmin)
