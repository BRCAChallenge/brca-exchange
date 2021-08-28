from django.core.management.base import BaseCommand, CommandError
from django.db import transaction
from data.utilities import update_findlay_functional_assays, update_materialized_view

class Command(BaseCommand):
    help = 'Updates findlay functional assays retroactively'

    @transaction.atomic
    def handle(self, *args, **options):
        update_findlay_functional_assays()
        update_materialized_view()
        print("Done!")
