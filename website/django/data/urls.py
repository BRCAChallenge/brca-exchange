from django.conf.urls import url

from . import views

urlpatterns = [

    url(r'^$', views.index, name="index"),
    url(r'^suggestions/$', views.autocomplete),
    url(r'^ga4gh/variants/search$', views.search_variants, name='search_variants'),
    url(r'^ga4gh/variants/(?P<variant_id>.+)$', views.get_variant, name = 'get_variant'),
    url(r'^ga4gh/variantsets/search', views.search_variant_sets, name='search_variant_sets'),
    url(r'^ga4gh/variantsets/(?P<variantSetId>.+)$', views.get_variant_set, name='get_variant_set'),
    url(r'^ga4gh/datasets/search', views.search_datasets, name='search_datasets'),
    url(r'^ga4gh/variantsets', views.varsetId_empty_catcher, name='empty_request_catcher'),
    url(r'^ga4gh/variants', views.empty_varId_catcher, name='empty_varId_catcher')

]
