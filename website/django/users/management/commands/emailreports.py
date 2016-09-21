import json
import os
import random
import md5
import requests
from urllib2 import HTTPError
from django.core.mail import EmailMultiAlternatives
from django.template import Context
from django.template.loader import get_template

from brca import settings, site_settings
from users.models import MyUser
from users.admin import UnapprovedUser

from django.core.management.base import BaseCommand, CommandError

class Command(BaseCommand):
    help = 'Send email notification to admins of users awaiting approval; intended for cron job'

    def handle(self, *args, **options):
        # Check for users awaiting approval
        users_awaiting_approval = MyUser.objects.filter(is_active = True, is_approved = False)
        user_count = users_awaiting_approval.count()

        if user_count == 0:
            self.stdout.write(self.style.SUCCESS("No users awaiting approval."))
        else:
            self.stdout.write("There are %i users awaiting approval." % user_count)

            recipients = MyUser.objects.filter(admin_notifications = True).filter(is_admin = True).values_list('email', flat=True)
            plaintext_email = get_template(os.path.join(settings.BASE_DIR, 'users', 'templates', 'users_awaiting_approval_email.txt'))
            html_email = get_template(os.path.join(settings.BASE_DIR, 'users', 'templates', 'users_awaiting_approval_email.html'))
            # get backend url somehow?
            url = "http://brcaexchange.org/backend/admin/users/unapproveduser/"
            d = Context({'user_count': user_count, 'url': url})

            subject, from_email = 'admin notification test', 'noreply@brcaexchange.org'
            text_content = plaintext_email.render(d)
            html_content = html_email.render(d)
            msg = EmailMultiAlternatives(subject, text_content, from_email, recipients)
            msg.attach_alternative(html_content, "text/html")
            msg.send()

            self.stdout.write(self.style.SUCCESS("email sent"))
