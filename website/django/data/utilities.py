from django.db import connection
from models import DataRelease, CurrentVariant, Variant, MupitStructure
import requests
import json
import sys
import time


def update_autocomplete_words():
        # Drop words table and recreate with latest data
        with connection.cursor() as cursor:
            cursor.execute("""
                DROP TABLE IF EXISTS words;
                CREATE TABLE words AS SELECT DISTINCT left(word, 300) as word, release_id FROM (
                SELECT regexp_split_to_table(lower("Genomic_Coordinate_hg38"), '[\s|''"]') as word, "Data_Release_id" as release_id from variant UNION
                SELECT regexp_split_to_table(lower("Genomic_Coordinate_hg37"), '[\s|''"]') as word, "Data_Release_id" as release_id from variant UNION
                SELECT regexp_split_to_table(lower("Genomic_Coordinate_hg36"), '[\s|''"]') as word, "Data_Release_id" as release_id from variant UNION
                SELECT regexp_split_to_table(lower("Clinical_significance_ENIGMA"), '[\s|''"]') as word, "Data_Release_id" as release_id from variant UNION
                SELECT regexp_split_to_table(lower("Gene_Symbol"), '[\s|''"]') as word, "Data_Release_id" as release_id from variant UNION
                SELECT regexp_split_to_table(lower("Reference_Sequence"), '[\s|''"]') as word, "Data_Release_id" as release_id from variant UNION
                SELECT regexp_split_to_table(lower("HGVS_cDNA"), '[\s|:''"]')  as word, "Data_Release_id" as release_id from variant UNION
                SELECT regexp_split_to_table(lower("BIC_Nomenclature"), '[\s|''"]') as word, "Data_Release_id" as release_id from variant UNION
                SELECT regexp_split_to_table(lower("HGVS_Protein"), '[\s|''"]') as word, "Data_Release_id" as release_id from variant
                )
                AS combined_words;

                CREATE INDEX words_idx ON words(word text_pattern_ops);
            """)


def set_release_name_defaults(apps, schema_editor):
    count = 1
    releases = DataRelease.objects.all().order_by('date')
    for release in releases:
        release.name = count
        release.save()
        count += 1


def send_mupit_request(query_url, params, retries=5):
    if retries <= 0:
        print "Request failed 5 times, exiting."
        sys.exit(1)
    try:
        r = requests.post(query_url, data=params)
        return r
    except requests.exceptions.RequestException as e:
        print e
        time.sleep(10)
        retries -= 1
        r = send_mupit_request(query_url, params, retries)
    return r


def update_mupit_structure_for_existing_variants(apps, schema_editor):
    mupit_structures = {ms['name']: ms['id'] for ms in MupitStructure.objects.values()}
    cvs = CurrentVariant.objects.all()
    for cv in cvs:
        variant = Variant.objects.get(id=getattr(cv, 'id'))
        chrom = "chr" + getattr(variant, 'Chr')
        pos = int(getattr(variant, 'Pos'))
        if (pos >= 32356427 and pos <= 32396972) or (pos >= 43045692 and pos <= 43125184):
            main_url = 'http://staging.cravat.us/MuPIT_Interactive'
            brca_structures = ['1t15','1jm7','4igk','fENSP00000380152_7']
            query_url = main_url+'/rest/showstructure/query'
            params = {
                     'search_textarea':'%s %s'%(chrom, pos),
                     'search_gene':'',
                     'search_structure':'',
                     'search_protein':'',
                     'search_upload_file':'',
                     }

            r = send_mupit_request(query_url, params)
            d = json.loads(r.text)
            structures = d['structures']
            main_struct = None;
            min_pref_level = sys.maxsize # max size integer
            max_no_query_pos = -1
            max_no_res = -1
            pref_level_hit = False
            for structure_id in structures:
                structure = structures[structure_id]
                no_query_pos = len(structure['gmtoseqres'])
                no_res = structure['nores']
                pref_level = structure['prefLevel']
                if (pref_level >= 1) and (pref_level < min_pref_level):
                    pref_level_hit = True
                    min_pref_level = pref_level
                    max_no_query_pos = no_query_pos
                    max_no_res = no_res
                    main_struct = structure_id
                elif (pref_level == min_pref_level) or (pref_level == 0 and not(pref_level_hit)):
                    if no_query_pos > max_no_query_pos:
                        max_no_query_pos = no_query_pos
                        max_no_res = no_res
                        main_struct = structure_id
                    elif no_query_pos == max_no_query_pos:
                        if no_res > max_no_res:
                            max_no_res = no_res
                            main_struct = structure_id
            if main_struct in brca_structures:
                mupit_structure_id = mupit_structures[main_struct]
                setattr(variant, "Mupit_Structure_id", mupit_structure_id)
                variant.save()
            time.sleep(0.1)
