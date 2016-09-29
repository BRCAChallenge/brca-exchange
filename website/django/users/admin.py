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

class ProxyUser(MyUser):
    class Meta:
        proxy = True
        verbose_name = "user"
        verbose_name_plural = "All Users"

class MyUserAdmin(admin.ModelAdmin):

    list_display = ('email', 'is_admin', 'firstName', 'lastName', 'is_approved', 'is_active', image)
    exclude = ('password',)
    ordering = ('is_approved',)
    search_fields = ('firstName', 'lastName', 'email')
    formfield_overrides = {
        models.TextField: {'widget': TextInput(attrs={'size':'20'})}
    }

    actions = ('approve',)

    def approve(self, request, queryset):
        queryset.update(is_approved=True)

class InactiveUser(MyUser):
    class Meta:
        proxy = True
        verbose_name = "inactive user"
        verbose_name_plural = "Users Awaiting Email Confirmation"

class InactiveUserAdmin(MyUserAdmin):
    def get_queryset(self, request):
        return self.model.objects.filter(is_active=False)

class UnapprovedUser(MyUser):
    class Meta:
        proxy = True
        verbose_name = "unapproved user"
        verbose_name_plural = "Users Awaiting Approval"

class UnapprovedUserAdmin(MyUserAdmin):
    def get_queryset(self, request):
        return self.model.objects.filter(is_approved=False, is_active=True)

admin.site.register(ProxyUser, MyUserAdmin)
admin.site.register(InactiveUser, InactiveUserAdmin)
admin.site.register(UnapprovedUser, UnapprovedUserAdmin)

