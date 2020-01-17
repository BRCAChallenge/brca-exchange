from django.db import connection
from .models import DataRelease, CurrentVariant, Variant, MupitStructure
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


def update_materialized_view():
    with connection.cursor() as cursor:
        cursor.execute(
            """
            DROP MATERIALIZED VIEW IF EXISTS currentvariant;
                CREATE MATERIALIZED VIEW currentvariant AS (
                    SELECT * FROM "variant" WHERE (
                        "id" IN ( SELECT DISTINCT ON ("Genomic_Coordinate_hg38") "id" FROM "variant" ORDER BY "Genomic_Coordinate_hg38" ASC, "Data_Release_id" DESC )
                    )
                );
            """
        )


def set_release_name_defaults(apps, schema_editor):
    count = 1
    releases = DataRelease.objects.all().order_by('date')
    for release in releases:
        release.name = count
        release.save()
        count += 1


def send_mupit_request(query_url, params):
    MAX_TRIES = 5
    tries = 0
    resp = None
    while tries < MAX_TRIES:
        try:
            resp = requests.post(query_url, data=params)
            return resp
        except requests.exceptions.RequestException as e:
            print(e)
            time.sleep(10)
            tries += 1
            continue
        break
    print("Request failed 5 times, exiting.")
    sys.exit(1)


def is_point_substitution(ref, alt):
    bases = ['a', 'c', 't', 'g']
    ref = ref.lower()
    alt = alt.lower()
    if ref in bases and alt in bases and len(ref) == 1 and len(alt) == 1:
        return True
    return False


def is_relevant_position(pos):
    # If positions requires updates, also update isRelevantPosition in pipeline/data_merging/getMupitStructure.py
    return (pos >= 32356427 and pos <= 32396972) or (pos >= 43045692 and pos <= 43125184)


def has_relevant_protein_change(hgvs_protein):
    # NOTE: if this requires updating, also update hasRelevantProteinChange in pipeline/data_merging/getMupitStructure.py
    # ignore variants with no amino acid change or variants with an unknown amino acid change
    if "?" in hgvs_protein or "=" in hgvs_protein or "*" in hgvs_protein:
        return False

    # remove all numbers and brackets from protein
    stripped_hgvs = ''.join([i for i in hgvs_protein if not (i.isdigit() or i in ['(', ')'])])
    hgvs_last_three_chars = stripped_hgvs[-3:]

    # ignore variants with an introduced stop codon
    if hgvs_last_three_chars == "Ter":
        return False

    return True


def update_mupit_structure_for_existing_variants():
    mupit_structures = {ms['name']: ms['id'] for ms in list(MupitStructure.objects.values())}
    cvs = CurrentVariant.objects.all()
    for cv in cvs:
        variant = Variant.objects.get(id=getattr(cv, 'id'))
        chrom = "chr" + getattr(variant, 'Chr')
        pos = int(getattr(variant, 'Pos'))
        ref = getattr(variant, 'Ref')
        alt = getattr(variant, 'Alt')
        hgvs_protein = getattr(variant, 'HGVS_Protein')
        if is_point_substitution(ref, alt) and is_relevant_position(pos) and has_relevant_protein_change(hgvs_protein):
            '''
            When you submit a position, the mupit server responds with a bunch of structures that the position maps to.
            Which structure to display first is determined client side.
            It chooses which structure to display based on three tests.
            If a higher test results in a tie, the next test is used to break the tie.
            1) pref_level. Select the lowest pref_level, but ignore pref_levels of zero
            2) no_query_pos. Select the structure that has the most hits from the submitted positions.
            Mupit supports querying with multiple positions at once, though brca exchange does not do this.
            3) no_res. Select the structure with the most residues (amino acids)
            '''
            main_url = 'http://mupit.icm.jhu.edu/MuPIT_Interactive'
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


class Benchmark(object):
    def __init__(self,name):
        self.name = name

    def __enter__(self):
        self.start = time.time()

    def __exit__(self,ty,val,tb):
        end = time.time()
        print(("%s : %0.3f seconds" % (self.name, end-self.start)))
        return False
