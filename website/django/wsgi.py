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

site.addsitedir(os.path.join(appdir, '../virtualenv/lib/python2.7/site-packages'))
from django.core.wsgi import get_wsgi_application

sys.path.append(appdir)

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "brca.settings")

application = get_wsgi_application()
