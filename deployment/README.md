# azure instances

Two VMs are provisioned: brcaexchange.cloudapp.net for production, and
brcaexchange-dev.cloudapp.net for development.

# virtual hosts

On each VM, two apache virtual hosts are configured: one for beta, and one for production.
The beta vhost is on the azure domain, (e.g. brcaexchange.cloudapp.net). The production vhost is
on the public domain (e.g. brcaexchange.org).

No public domain has been published for brcaexchange-dev.cloudapp.net, but you
may add the following line to your /etc/hosts file if you need access to the
'production' virutal host:

40.78.99.255 brcaexchange-dev.cloudapp.net

# apache directories
## beta vhost

| --- | --- |
| /var/www/html/beta | static assets |
| /var/www/backend/beta/virtualenv | the python virtualenv |
| /var/www/backend/beta/django | application django source |

## production vhost

| --- | --- |
| /var/www/html/production | static assets |
| /var/www/backend/production/virtualenv | the python virtualenv |
| /var/www/backend/production/django | application django source |

# apache configuration

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

## javascript

An apache server-side include (SSI) is used to inject the configuration into
index.html. The SSI 'include' statement is in website/page.template, which is
the template used to generate index.html. The apache config includes directives
to enable SSI on html files in the static asset directory.

```
    <Directory /var/www/html/beta>
	Options +Includes
	AddOutputFilter INCLUDES .html
    </Directory>
```
# Continuous integration

CI is implemented with circle-ci (CCI). CCI will build any commits on 'master' branch,
and deploy them on the beta vhost of brcaexchange-dev. The deploy script, deploy-dev,
will rsync the build, and copy site settings into place.
