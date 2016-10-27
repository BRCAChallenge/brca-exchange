from django.conf.urls import url

from . import views

urlpatterns = [

    url(r'^$', views.index, name="index"),
    url(r'^releases', views.releases, name="releases"),
    url(r'^variant/$', views.variant, name="variant"),
    url(r'^suggestions/$', views.autocomplete),
    url(r'^ga4gh/v0.6.0a7/variants/search$', views.search_variants, name='search_variants'),
    url(r'^ga4gh/v0.6.0a7/variants/(?P<variant_id>.+)$', views.get_variant, name = 'get_variant'),
    url(r'^ga4gh/v0.6.0a7/variantsets/search', views.search_variant_sets, name='search_variant_sets'),
    url(r'^ga4gh/v0.6.0a7/variantsets/(?P<variant_set_id>.+)$', views.get_variant_set, name='get_variant_set'),
    url(r'^ga4gh/v0.6.0a7/datasets/search', views.search_datasets, name='search_datasets'),
    url(r'^ga4gh/v0.6.0a7/datasets/(?P<dataset_id>.+)$', views.get_dataset, name="get_dataset" ),
    url(r'^ga4gh/v0.6.0a7/variantsets', views.empty_variantset_id_catcher, name='empty_variantset_id_catcher'),
    url(r'^ga4gh/v0.6.0a7/variants', views.empty_variant_id_catcher, name='empty_variant_id_catcher'),
    url(r'^ga4gh/v0.6.0a7/datasets', views.empty_dataset_catcher, name='empty_dataset_catcher')

]
