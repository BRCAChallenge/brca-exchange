from django.conf.urls import url

from . import views

urlpatterns = [

    url(r'^$', views.index, name="index"),
    url(r'^suggestions/$', views.autocomplete),
    url(r'^ga4gh/variants/search$', views.variant_search, name='variant_search'),
    url(r'^ga4gh/variants/(?P<variant_id>.+)$', views.get_var_by_id, name = 'get_var_by_id'),
    url(r'^ga4gh/variantsets/search', views.get_variantSet, name='get_variantSet'),
    url(r'^ga4gh/variantsets/(?P<variantSetId>.+)$', views.get_varset_by_id, name='get_varset_by_id'),
    url(r'^ga4gh/datasets/search', views.search_datasets, name='search_datasets'),
    url(r'^ga4gh/variantsets', views.varsetId_empty_catcher, name='empty_request_catcher'),
    url(r'^ga4gh/variants', views.empty_varId_catcher, name='empty_varId_catcher')

]
