from django.core.management.base import BaseCommand, CommandError
from django.db import transaction
from data.utilities import update_mupit_structure_for_existing_variants

class Command(BaseCommand):
    help = 'Updates mupit structures'

    @transaction.atomic
    def handle(self, *args, **options):
        update_mupit_structure_for_existing_variants()
        print "Done!"
