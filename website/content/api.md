# BRCA Exchange GA4GH API

<div style="display:inline-block;text-align:center">
    <a href="http://genomicsandhealth.org"><img src="ga4gh-logo-more.png" style="width:240px;"></a>
    <a href="https://genomics.soe.ucsc.edu"><img src="ucsc_logo.png" style="width:260px;"></a>
</div>

---

The GA4GH has defined a standard interface for exchanging genomic data over HTTP called an API, or Application Programming Interface.

Using the <a href="http://github.com/ga4gh/schemas">GA4GH schemas</a>, students at UCSC have developed an implementation of this interface that provides programmatic access to the BRCA Exchange variation data and annotations. 

This allows the well curated and expert reviewed variation data in the BRCA exchange to be made available to application developers in most common languages.

The BRCA exchange exposes the interfaces necessary for accessing variant data, which are a subset of the interfaces the GA4GH defines. For a detailed description of the interface view the generated <a href="/backend/static/swagger/index.html"><b>documentation</b></a>.

An example of accessing variants through GA4GH client in a more systematica way can be obtained through the <a href="https://github.com/achave11/brca-exchange/blob/master/website/django/python_notebook/brca-exchange.ipynb"><b>python notebook</b></a>.