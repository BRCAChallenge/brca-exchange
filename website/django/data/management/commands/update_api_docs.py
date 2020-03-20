import os, re, subprocess
import sys
from argparse import FileType

from django.core.management import BaseCommand
from django.template.loader import render_to_string


class Command(BaseCommand):
    help = 'Processes the API documentation template'

    def add_arguments(self, parser):
        parser.add_argument(
            'outfile', type=FileType('w'), default=sys.stdout, nargs='?',
            help='Where to put the processed API docs'
        )

    def handle(self, *args, **options):
        import data
        app_path = os.path.dirname(data.__file__)

        views_path = os.path.join(app_path, "views.py")

        # open django/data/views.py and extract all function names
        func_re = re.compile(r'^def ([A-Za-z][^(]+)\(')
        functions = dict(
            (func_re.match(line).groups()[0], idx+1)
            for idx, line in enumerate(open(views_path))
            if func_re.match(line)
        )

        # also get the current commit hash so we can embed it in the github URL
        # commit = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"])
        # alternatively, we can get the most recent commit to django/data/views.py...
        commit = subprocess.check_output(("git log -n 1 --pretty=format:%%h -- %s" % views_path).split(" "))

        options['outfile'].write(
            render_to_string('api_overview.md', {
                'views': functions,
                'commit': commit.decode("utf-8").strip()
            })
        )
