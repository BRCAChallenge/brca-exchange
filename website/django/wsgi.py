"""
WSGI config for brca project.

It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/1.9/howto/deployment/wsgi/
"""

import os
import sys
import site

appdir = os.path.dirname(__file__)

site.addsitedir(os.path.join(appdir, '../virtualenv3/lib/python3.6/site-packages'))
from django.core.wsgi import get_wsgi_application
from django.db.backends.signals import connection_created
from django.dispatch import receiver

sys.path.append(appdir)

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "brca.settings")

@receiver(connection_created)
def setup_postgres(connection, **kwargs):
    if connection.vendor != 'postgresql':
        return
    # Timeout statements after 45 seconds.
    with connection.cursor() as cursor:
        cursor.execute("""
            SET statement_timeout TO 45000;
        """)


application = get_wsgi_application()
