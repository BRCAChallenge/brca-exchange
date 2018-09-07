from django.core.management.base import BaseCommand, CommandError
# from django.conf import settings
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
        variant_papers = literature_results['variants']
        papers = literature_results['papers']
        paper_objects = {}
        for pmid, paper in papers.iteritems():
            query = Paper.objects.filter(pmid=pmid)
            if query.count() > 0:
                # we already have this paper in the database
                paper_objects[pmid] = query[0]
            else:
                p = Paper(title=paper['title'], authors=paper['authors'], journal=paper['journal'], \
                        keywords=paper['keywords'], abstract=paper['abstract'], year=paper['year'], \
                        pmid=paper['pmid'])
                p.save()
                paper_objects[pmid] = p
        for variant_name, variant_papers in variant_papers.iteritems():
            for variant_paper in variant_papers:
                pmid = variant_paper['varId']
                if pmid in paper_objects:
                    paper = paper_objects[pmid]
                    mentions = variant_paper['geneSnippets']
                    if mentions == None:
                        mentions = ''
                    vp = VariantPaper(variant_hg38=variant_name, paper=paper, mentions=mentions)
                    vp.save()
