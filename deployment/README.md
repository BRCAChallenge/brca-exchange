# Release process

See [this page](RELEASES.md).

# azure instances

Two VMs are provisioned: brcaexchange.cloudapp.net for production, and
brcaexchange-dev.cloudapp.net for development.

# virtual hosts

On each VM, two apache virtual hosts are configured: one for beta, and one for production.
The beta vhost is on the azure domain, (e.g. brcaexchange.cloudapp.net). The production vhost is
on the public domain (e.g. brcaexchange.org).

| VM | vhost | domain |
| --- | --- | --- |
| brcaexchange |  beta | brcaexchange.cloudapp.net |
| brcaexchange |  production | brcaexchange.org |
| brcaexchange-dev |  beta | brcaexchange-dev.cloudapp.net |
| brcaexchange-dev |  production | |


No public domain has been published for brcaexchange-dev.cloudapp.net, but you
may add the following line to your /etc/hosts file if you need access to the
'production' virtual host on brcaexchange-dev.cloudapp.net:

```
40.78.99.255 brcaexchange-dev.org
```

# apache directories
## beta vhost

| file | description |
| --- | --- |
| /var/www/html/beta | static assets |
| /var/www/backend/beta/virtualenv | the python virtualenv |
| /var/www/backend/beta/django | application django source |

## production vhost

| file | description |
| --- | --- |
| /var/www/html/production | static assets |
| /var/www/backend/production/virtualenv | the python virtualenv |
| /var/www/backend/production/django | application django source |

# apache configuration

| file | description |
| --- | --- |
| sites-available/000-default | configuration for the production site |
| sites-available/001-beta | configuration for the beta site |

# site settings

django and javascript settings are required for each vhost, to set keys
and database details. These are stored in ```~brca/site_settings```, and
must be copied into the application directories during deployment. E.g.

```
cp ~brca/site_settings/config.beta.js /var/www/html/beta/config.js
cp ~brca/site_settings/site_settings.beta.py /var/www/backend/beta/django/brca/site_settings.py
```

This is done by the ```deploy-dev``` and ```deploy-production``` scripts. ```deploy-dev``` will
be run automatically by the CI (see below). ```deploy-production``` is run manually at the end
of QA.

## javascript

An apache server-side include (SSI) is used to inject the configuration into
```index.html```. The SSI 'include' statement is in ```website/page.template```, which is
the template used to generate ```index.html```. The apache config includes directives
to enable SSI on html files in the static asset directory.

```
    <Directory /var/www/html/beta>
	Options +Includes
	AddOutputFilter INCLUDES .html
    </Directory>
```

# Continuous integration

CI is implemented with circle-ci 2.0 (CCI). CCI will build any commits on 'master' branch,
and deploy them on the beta vhost of brcaexchange-dev. The deploy script, deploy-dev,
will rsync the build, and copy site settings into place.

Any commits that are tagged with a version number, e.g. v1.2.3, will be built and
deployed to brcaexchange.cloudapp.net (beta). The ```release``` script will
generate and increment tags of this form automatically.

To allow CCI to push code to azure, an authorized ssh key is added to the ssh config
on the azure hosts. The ssh key can be managed in the CCI account settings.

You may log in to circle-ci with your github credentials.

## Running Circle CI Locally

Version 2.0 of CCI can easily be run locally as well using

```
circleci build
```

See the [local client documentation](https://circleci.com/docs/2.0/local-cli/) for more details, also regarding installation.

# Community Data

##### WARNING: This will delete and replace any relevant data in the beta database and `/var/www/backend/beta/django/uploads/` with data from production. Affected tables include all `auth` prefixed tables, all `users` prefixed tables, `django_admin_log`, and `django_content_type`.

To copy community data from the production server and database to beta, run `sh deployment/copy-community-data.sh`.
