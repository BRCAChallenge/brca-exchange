from django.core.management.base import BaseCommand, CommandError
from django.db import connection, transaction
from data.models import CurrentVariant, Variant, Paper, VariantPaper
from argparse import FileType
import json
from django.core.management import call_command



class Command(BaseCommand):
    help = 'Add literature search results to database'

    def add_arguments(self, parser):
        parser.add_argument('literature', type=FileType('r'), help='Literature search results to add, in JSON format')

    @transaction.atomic
    def handle(self, *args, **options):
        literature_results = json.load(options['literature'])
        variants_found_in_papers = literature_results['variants']
        papers = literature_results['papers']

        # Soft delete all existing records (they will be undeleted if they're
        # in the new data)
        existing_papers = Paper.objects.all()
        existing_papers.update(deleted=True)
        existing_variant_papers = VariantPaper.objects.all()
        existing_variant_papers.update(deleted=True)

        paper_objects = {}
        for pmid, paper in papers.iteritems():
            query = Paper.objects.filter(pmid=pmid)
            if query.count() > 0:
                # we already have this paper in the database
                paper_objects[pmid] = query[0]
                query.update(deleted=False)
            else:
                p = Paper(title=paper['title'], authors=paper['authors'], journal=paper['journal'], \
                        keywords=paper['keywords'], abstract=paper['abstract'], year=paper['year'], \
                        deleted=False, pmid=paper['pmid'])
                p.save()
                paper_objects[pmid] = p

        for variant_genomic_coordinate, variant in variants_found_in_papers.iteritems():
            for pmid, mentions in variant.iteritems():
                if pmid in paper_objects:
                    paper = paper_objects[pmid]
                    if mentions == None:
                        mentions = []
                    query = VariantPaper.objects.filter(variant_hg38=variant_genomic_coordinate, paper=paper)
                    if query.count() > 0:
                        # we already have this variantpaper
                        query.update(mentions=mentions, deleted=False)
                    else:
                        vp = VariantPaper(variant_hg38=variant_genomic_coordinate, paper=paper, mentions=mentions, deleted=False)
                        vp.save()
