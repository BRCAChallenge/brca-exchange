#!/bin/bash
set -o nounset
set -o errexit

# NOTE: This script copies community data from production to beta.

######################################################################################
#   WARNING: All `users` and `auth` prefixed tables, as well as `django_admin_log`   #
#   `and django_content_type` tables will be dropped and all data will be replaced.  #
######################################################################################

HOST=${HOST:-brcaexchange.cloudapp.net}
USER=brca

# Copy community relevant DB data from production to beta
ssh -l${USER} ${HOST} <<-ENDSSH
    set -o errexit
    sudo chown -R www-data:www-data /var/www/backend/beta/django/uploads
    . /var/www/backend/beta/virtualenv/bin/activate
    sudo rm -rf /var/www/backend/beta/django/uploads
    cp -r /var/www/backend/production/django/uploads /var/www/backend/beta/django/uploads
    sudo -u postgres pg_dump -d production.pg -F c -t "^users*" -t django_admin_log -t django_content_type -t "^auth*" -c -f /tmp/users_seq.dump
    sudo -u postgres pg_restore /tmp/users_seq.dump -c -v -1 -d storage.pg
    sudo apache2ctl restart
ENDSSH
