from django.urls import re_path

from . import views

urlpatterns = [

    re_path(r'^$', views.index, name="index"),
    re_path(r'^releases', views.releases, name="releases"),
    re_path(r'^variant/$', views.variant, name="variant"),
    re_path(r'^vrid/$', views.vrid, name="vr_id"),
    re_path(r'^variant/(?P<variant_id>.+)/reports$', views.variant_reports, name="variant_reports"),
    re_path(r'^variantcounts', views.variant_counts, name="variant_counts"),
    re_path(r'^variantpapers/$',views.variant_papers, name="variant_papers"),
    re_path(r'^variantreps/$', views.variantreps, name="variant_reps"),
    re_path(r'^sitemap.txt$', views.sitemap, name="sitemap"),
    re_path(r'^suggestions/$', views.autocomplete)

]
