"""
WSGI config for brca project.

It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/1.9/howto/deployment/wsgi/
"""

import os
import sys

from django.core.wsgi import get_wsgi_application

appdir = '/var/www/backend/beta/django'
sys.path.append(appdir)

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "brca.settings")

application = get_wsgi_application()
