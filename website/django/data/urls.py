from django.conf.urls import url

from . import views

urlpatterns = [

    url(r'^$', views.index, name="index"),
    url(r'^releases', views.releases, name="releases"),
    url(r'^variant/$', views.variant, name="variant"),
    url(r'^variant/(?P<variant_id>.+)/reports$', views.variant_reports, name="variant_reports"),
    url(r'^variantcounts', views.variant_counts, name="variant_counts"),
    url(r'^variantpapers/$',views.variant_papers, name="variant_papers"),
    url(r'^variantreps/$', views.variantreps, name="variant_reps"),
    url(r'^sitemap.txt$', views.sitemap, name="sitemap"),
    url(r'^suggestions/$', views.autocomplete)

]
