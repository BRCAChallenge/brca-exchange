# Internal API Documentation

This document describes BRCA Exchange's internal API, which is used to serve data to the web and mobile frontends. Notice that the API is subject to change at any time and without warning, but this document will be kept up-to-date as changes occur.

All paths are prefixed with `https://brcaexchange.org/backend/data`.

## Overview

Endpoints are listed in the following table.

| Path                                       | Name                                                                    | Code Reference |
| ------------------------------------------ | ----------------------------------------------------------------------- | -------------- |
| /                                          | [index](#user-content---search-index)                                   | [views.index](https://github.com/BRCAChallenge/brca-exchange/blob/444110d6/website/django/data/views.py#L258) |
| /releases                                  | [releases](#user-content-releases--release--list)                       | [views.releases](https://github.com/BRCAChallenge/brca-exchange/blob/444110d6/website/django/data/views.py#L30) |
| /variant/                                  | [variant](#user-content-variant-variant-details)                        | [views.variant](https://github.com/BRCAChallenge/brca-exchange/blob/444110d6/website/django/data/views.py#L94) |
| /vrid/                                     | [vrid](#user-content-vrid-variant-representation-vr-id)                 | [views.vrid](https://github.com/BRCAChallenge/brca-exchange/blob/444110d6/website/django/data/views.py#L) |
| /variant/(&lt;variant_id&gt;)/reports      | [variant reports](#user-content-variantvariant_idreports-reports)       | [views.variant_reports](https://github.com/BRCAChallenge/brca-exchange/blob/444110d6/website/django/data/views.py#L151) |
| /variantcounts                             | [variant counts](#user-content-variantcounts-variant-counts)            | [views.variant_counts](https://github.com/BRCAChallenge/brca-exchange/blob/444110d6/website/django/data/views.py#L50) |
| /variantpapers                             | [variant papers](#user-content-variantpapers--literature-results)       | [views.variant_papers](https://github.com/BRCAChallenge/brca-exchange/blob/444110d6/website/django/data/views.py#L189) |
| /variantreps/                              | [variant reps](#user-content-variantreps--variant-representation-table) | [views.variantreps](https://github.com/BRCAChallenge/brca-exchange/blob/444110d6/website/django/data/views.py#L111) |
| /sitemap.txt                               | [sitemap](#user-content-sitemaptxt-sitemap-metadata)                    | [views.sitemap](https://github.com/BRCAChallenge/brca-exchange/blob/444110d6/website/django/data/views.py#L131) |
| /suggestions/                              | [suggestions](#user-content-suggestions-suggestions)                    | [views.autocomplete](https://github.com/BRCAChallenge/brca-exchange/blob/444110d6/website/django/data/views.py#L582) |

*(Pay special attention to each path's trailing slash (or lack thereof). If the slash is included, it must be supplied, and omitted otherwise.)*


## Endpoint Descriptions

<!-- -------------------------------------------------------------------- -->
### /:  Search Index
<!-- -------------------------------------------------------------------- -->

Search endpoint; returns variants based on ordering and filtering criteria.

#### Request

| Parameter      | Description                                                                              | Default |
|----------------|------------------------------------------------------------------------------------------|---------|
| order_by       | column by which to order results                                                         | None    |
| direction      | 'descending' for descending ordering, or anything else for ascending                     | None    |
| page_size      | Number of results in the page; 0 to disable paging, which returns all variants           | 0       |
| page_num       | Index of the current page                                                                | 0       |
| search_term    | String to search for in column values; unused if None                                    | None    |
| format         | 'csv', 'tsv', or 'json' (server error if None)                                           | None*   |
| include[]      | Sources to include in results (OR'd); no results if None                                 | None*   |
| exclude[]      | Sources to exclude from results (OR'd)                                                   | None    |
| filter[]       | Columns to filter according to the values in `filterValues`                              | None    |
| filterValue[]  | Values which must occur in the columns specified in `filter`                             | None    |
| column[]       | Unused?                                                                                  | None    |
| release        | ID of the release from which to return results; any release if None                      | None    |
| change_types[] | Only used if `release` is specified; returns variants with the given change types (OR'd) | None    |
| show_deleted   | whether to include deleted variants in the results; True if any value is specified       | False   |

*( * required; the query won't complete successfully unless a value is specified for these parameters. )*

Parameters that end with "[]" are arrays; you can specify the same parameter repeatedly to add entries to the list. The order of elements
matters in the case of `filter`/`filterValue`; corresponding elements are zipped together to produce per-column filtering predicates. For example,

`?filter=Protein_Change&filter=Pathogenicity_expert&filterValue=V600E&filterValue=Pathogenic`

filters for variants with the `Protein_Change` column containing the string "V600E" and the `Pathogenicity_expert` column containing the string "Pathogenic".

Special notes about parameters:
- `format`: Required, produces a server error if unspecified. See Response below for information on how this parameter affects the response.
- `search_term`: If specified, this value must occur in at least one column for the variant to be returned. The logic is a little complicated; most, but not all fields, are searched, and most are searched **without** case sensitivity. The code is in [apply_search()](https://github.com/BRCAChallenge/brca-exchange/blob/master/website/django/data/views.py#L375). Synonym matching is also handled in that routine, affecting the value of the `synonyms` field in the JSON response.
- `filter[]`/`filterValue[]`: Pairs of a column name and the (case-sensitive) value which must occur in that column.
- `include[]`: names of columns which must be the literal value True in the result set, restricting the use of this parameter to boolean-typed columns. In the web UI, the choices for this parameter are populated from the "source" columns [defined in the Variant model](https://github.com/BRCAChallenge/brca-exchange/blob/master/website/django/data/models.py#L64). Possible values for this parameter are:
  - Variant_in_ENIGMA
  - Variant_in_ClinVar
  - Variant_in_1000_Genomes
  - Variant_in_ExAC
  - Variant_in_LOVD
  - Variant_in_BIC
  - Variant_in_ESP
  - Variant_in_exLOVD
  - Variant_in_Findlay_BRCA1_Ring_Function_Scores
  - Variant_in_GnomAD
- `exclude[]`: opposite of `include[]`; if the variant has a literal value of False for any of these columns, the variant is excluded from the results. Possible values are the same as for `include[].`
- `change_types[]`: if `release` is specified, variants with the types of changes specified in this array parameter will be included in the results. Possible for this parameter are:
  - new
  - deleted
  - added_classification
  - changed_classification
  - added_information
  - changed_information
  - none


#### Response

If `format` is "json", data is returned in this form:

```
{
  count: <num_total:int>
  data: [<:Variant>, ...]
  deletedCount: <num_deleted:int>
  releaseName: <release_name:string>|null
  synonyms: <num_synonyms:int>
}
```

- `count`: the total number of returned variants, disregarding paging.
- `deletedCount`: the number of results that would have been returned had they not been deleted. This field is null if `show_deleted` or `release` is specified.
- `releaseName` if `release` is specified, the name of returned release, and null otherwise.
- `data`: list of returned results, starting from page `page_num` and up to `page_size` in length.
- `synonyms`: If `search_term` is specified, the number of variants matched as synonyms of the query (i.e., the value matched something in the `Synonyms` or  `BIC_Nomenclature` columns).

"Variant" entries in `data` are instances of the [Variant](https://github.com/BRCAChallenge/brca-exchange/blob/master/website/django/data/models.py#L62) model. The `id` field indicates the ID of this version of the variant, which is used elsewhere in the API to query for more information about this variant. Each version of a variant receives a unique ID. An extra field, `Data_Release_id`, appears for each result, indicating the release which contains the most recent data for this variant. This ID corresponds to instances of the [DataRelease](https://github.com/BRCAChallenge/brca-exchange/blob/master/website/django/data/models.py#L7) model.

If `format` is "csv" or "tsv", `page_size` and `page_num` are ignored and all variants matching the query are returned as comma-separated or tab-separated values, respectively. Rows are delimited by newlines.


<!-- -------------------------------------------------------------------- -->
### /releases:  Release  List
<!-- -------------------------------------------------------------------- -->

Returns metadata about each release of the BRCA Exchange database.

#### Request

| Parameter      | Description                                                      | Default |
|----------------|------------------------------------------------------------------|---------|
| release_id     | ID of the release to retrieve; all releases are returned if None | None    |


#### Response

Data is returned in this form:

```
{
  latest: <release_id>
  releases: [<:Release>, ...]
}
```

- `latest`: The ID of the latest release; returned regardless of `release_id`
- `releases`: a list of [DataRelease](https://github.com/BRCAChallenge/brca-exchange/blob/master/website/django/data/models.py#L7) instances. If `release_id` is specified, this list will contain at most one release matching the ID, or no elements if no release with that ID is found.


<!-- -------------------------------------------------------------------- -->
### /variant: Variant Details
<!-- -------------------------------------------------------------------- -->

Returns detailed information about a single variant.

#### Request

| Parameter      | Description                                                      | Default |
|----------------|------------------------------------------------------------------|---------|
| variant_id     | ID of the variant to retrieve; server error if None              | None*   |


#### Response

```
{
  data: [<:Variant>, ...]
}
```

The array in `data` contains each version of the variant in reverse-chronological order (i.e., the most recent version is first). The first version
is the release in which this variant was first included in the database, and entries only exist when the variant changed between releases.

In addition to the fields in the Variant model, there are two extra fields in each variant: `Diff` and `Data_Release` (which replaces the `Data_Release_id` column from the search results endpoint). `Diff` contains the set of changes in the variant from the previous version, and is null for the first version. `Data_Release` is an object containing detailed information about the release for this version of the variant; it's the same data as in the `/releases` endpoint.

Diff entries have the following form:

```
{
  "field": <column changed:string>,
  "removed": <previous value of the column>,
  "added": <new value of the column>,
  "field_type": "individual"|"list"
}
```

If `field_type` is `individual`, `removed` and `added` contain scalar types. If `field_type` is `list`, then those fields will contain lists (or null). Throughout the system, '' (the empty string), '-', 'None', and null are used (somewhat) interchangeably as indicators of an empty field. You'll often see one of these values for `removed` if a value for that field was introduced in the current version.


<!-- -------------------------------------------------------------------- -->
### /vrid: Variant Representation (VR) ID
<!-- -------------------------------------------------------------------- -->

Given a VR ID, returns data about the associated variant. If the VR ID doesn't match any variant, the server
returns an HTTP 500 error.

#### Request

| Parameter      | Description                                                      | Default |
|----------------|------------------------------------------------------------------|---------|
| vr_id          | The VR hash for the desired variant                              | None*   |


#### Response

```
{
  data: [<:Variant>, ...]
}
```

The data returned in the Variant array is the same as returned from the `/variant` endpoint.


<!-- -------------------------------------------------------------------- -->
### /variant(&lt;variant_id&gt;)/reports: Reports
<!-- -------------------------------------------------------------------- -->

Returns reports associated with the variant with ID `variant_id`. Reports are the main source of information about the clinical impact of variants.

#### Request

This endpoint contains a path parameter, `variant_id`, which must be an integer identifying a valid variant ID. If no matching variant is found, the `data` list in the response will be empty. For ClinVar and  LOVD, historical reports (i.e., reports for releases before the release of the given variant ID) are also returned; these historical reports are matched by the key used by that source. ClinVar uses the `SCV_ClinVar` field and LOVD uses `Submission_ID_LOVD` to match historical reports.

#### Response

```
{
  "data": [<:Report>, ...]
}
```

The response is a list of [Report](https://github.com/BRCAChallenge/brca-exchange/blob/master/website/django/data/models.py#L350) instances, which are nearly identical to Variant instances, except for a few extra columns indicating the source, authors, date on which the report was generated, and other relevant metadata. As of this writing, the following sources exist in the database: ClinVar, exLOVD, ExAC, ENIGMA, Findlay_BRCA1_Ring_Function_Scores, LOVD, 1000_Genomes, ESP, BIC.


<!-- -------------------------------------------------------------------- -->
### /variantcounts: Variant Counts
<!-- -------------------------------------------------------------------- -->

Returns summary information about the contents of the database.

#### Request

This endpoint takes no parameters.

#### Response

The response is pretty self-explanatory. Sample results as of this writing (October 1st, 2019):

```json
{
  "brca1": {
    "pathogenic": 2189,
    "benign": 605,
    "total": 12628,
    "likelyPathogenic": 0,
    "likelyBenign": 441
  },
  "brca2": {
    "pathogenic": 2637,
    "benign": 639,
    "total": 13747,
    "likelyPathogenic": 0,
    "likelyBenign": 738
  },
  "enigma": 7256,
  "enigmaLikelyPathogenic": 0,
  "total": 26375,
  "enigmaBenign": 1244,
  "enigmaPathogenic": 4826,
  "enigmaLikelyBenign": 1179
}
```


<!-- -------------------------------------------------------------------- -->
### /variantpapers:  Literature Results
<!-- -------------------------------------------------------------------- -->

Returns literature associated with this variant (i.e., matching the variant's name, description, HGVS string, synonyms, etc.)

#### Request

| Parameter      | Description                                                            | Default |
|----------------|------------------------------------------------------------------------|---------|
| variant_id     | ID of the variant for which to retrieve results; server error if None  | None *  |


#### Response

Data is returned in this form:

```
{
  data: [<:VariantPaper + Paper>, ...]
}
```

- `data`: a list of [VariantPaper](https://github.com/BRCAChallenge/brca-exchange/blob/master/website/django/data/models.py#L632) instances for the given variant. The  VariantPaper object is extended by its corresponding [Paper](https://github.com/BRCAChallenge/brca-exchange/blob/master/website/django/data/models.py#L619) instance. (The same paper may mention multiple variants, thus papers and variant mentions within that paper are stored separately.)

Example: https://brcaexchange.org/backend/data/variantpapers/?variant_id=219222


<!-- -------------------------------------------------------------------- -->
### /variantreps/:  Variant Representation Table
<!-- -------------------------------------------------------------------- -->

For all variants in the current release, returns their VR representation (if available).

#### Request

This endpoint takes no parameters.

#### Response

Data is returned in this form:

```
{
  count: <num of results>
  data: [<:VariantRepresentation>, ...]
}
```

- `data`: a list of [VariantRepresentation](https://github.com/BRCAChallenge/brca-exchange/blob/master/website/django/data/models.py#L1019) instances.



<!-- -------------------------------------------------------------------- -->
### /sitemap.txt: Sitemap Metadata
<!-- -------------------------------------------------------------------- -->

Returns newline-delimited strings of crawlable URLs for the site. This endpoint is mostly so webcrawlers that require a sitemap (e.g., Google) can crawl the site effectively.

#### Request

This endpoint takes no parameters.

#### Response

Newline-delimited lines, e.g.:

```
https://brcaexchange.org/
https://brcaexchange.org/factsheet
https://brcaexchange.org/help
https://brcaexchange.org/community
https://brcaexchange.org/variants
https://brcaexchange.org/about/thisSite
https://brcaexchange.org/releases
https://brcaexchange.org/variant/146096
https://brcaexchange.org/variant/202542
https://brcaexchange.org/variant/215771
https://brcaexchange.org/variant/239814
https://brcaexchange.org/variant/199359
https://brcaexchange.org/variant/222196
```

The list begins with static URLs of pages on the site, then consists of URLs for the most recent version of each variant in the database (i.e., the contents of the `currentversion` view).


<!-- -------------------------------------------------------------------- -->
### /suggestions/: Suggestions
<!-- -------------------------------------------------------------------- -->

Quickly returns variants that partially match the given query term. This endpoint is used by the search box in the web UI to provide suggestions of completed values for partial queries before the user submits the query.

#### Request


https://brcaexchange.org/backend/data/suggestions/?term=c.17

| Parameter      | Description                                                                          | Default |
|----------------|--------------------------------------------------------------------------------------|---------|
| term           | the partial query so far; empty result if None                                       | None*   |
| release        | the release ID which to search for completions; uses the most current if unspecified | None    |
| limit          | maximum number of suggestions to return                                              | 10      |

#### Response

Suggestions are returned in this form:

```
{
  "suggestions": [[<completion>], ...]
}
```

`suggestions` contains each suggested completion of the partial query. Each suggestion is a string that would complete the partial query, wrapped in a list. The total number of elements will not exceed `limit`.

Suggestions are extracted as lowercase versions of distinct values that occur in the following columns in Variant:
- Genomic_Coordinate_hg38
- Genomic_Coordinate_hg37
- Clinical_significance_ENIGMA
- Gene_Symbol
- Reference_Sequence
- HGVS_cDNA
- BIC_Nomenclature
- HGVS_Protein

The suggestions are retrieved from the `words` table; the full logic for generating this table can be found here: [update_autocomplete_words()](https://github.com/BRCAChallenge/brca-exchange/blob/master/website/django/data/utilities.py#L9). This table is populated when a new release is loaded into the webserver's database via the `addrelease` management command.
