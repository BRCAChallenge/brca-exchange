"""
Backfill DataRelease.variants_added/variants_modified/variants_deleted from
production's live, publicly-computed numbers.

The finalized schema's Variant is current-state-only (no per-release history,
no Change_Type), so these can no longer be computed locally the way production
still computes them (see releases() on the master branch). Production's public
API is the only remaining source of truth for the historical values, so this
pulls them from there rather than reconstructing history that no longer exists
in this schema.

variants_classified is intentionally left at its default (0) — not backfilled.

Usage:
    python manage.py backfill_release_stats
    python manage.py backfill_release_stats --source https://brcaexchange.org/backend
"""

import requests

from django.core.management.base import BaseCommand
from data.models import DataRelease


class Command(BaseCommand):
    help = "Backfill variants_added/variants_modified/variants_deleted on DataRelease from production's live API"

    def add_arguments(self, parser):
        parser.add_argument(
            '--source', default='https://brcaexchange.org/backend',
            help='Base backend URL to pull release stats from (default: production)',
        )

    def handle(self, *args, **options):
        source = options['source'].rstrip('/')
        resp = requests.get(f'{source}/data/releases', timeout=60)
        resp.raise_for_status()
        stats_by_id = {r['id']: r for r in resp.json()['releases']}

        updated, missing = 0, []
        for release in DataRelease.objects.all():
            stats = stats_by_id.get(release.id)
            if stats is None:
                missing.append(release.id)
                continue
            release.variants_added = stats['variants_added']
            release.variants_modified = stats['variants_modified']
            release.variants_deleted = stats['variants_deleted']
            release.save(update_fields=['variants_added', 'variants_modified', 'variants_deleted'])
            updated += 1

        self.stdout.write(self.style.SUCCESS(f'Updated {updated} release(s) from {source}'))
        if missing:
            self.stdout.write(self.style.WARNING(
                f'{len(missing)} local release id(s) not found in source, left at 0: {missing}'
            ))
