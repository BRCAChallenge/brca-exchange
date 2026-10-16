# website for BRCA challenge
To contribute:
Fork the repository, make the changes, and submit pull request.

## Build the frontend
The build is based on npm and webpack.
* **Ensure that git and node are installed**
   * Node >=22.21.1 is required (see `engines` in `package.json`). Use [nvm](https://github.com/nvm-sh/nvm) if you need to manage multiple versions.
* **Start the frontend**
   * `cd website`
   * `npm install`
   * `npm start`
   * This starts webpack-dev-server on [http://localhost:8080/](http://localhost:8080/). Note that it does *not* point at a local backend by default — see "Point the frontend at a local backend" below.

## Build the server
The server runs on Django with Postgres, so install and set those up first.

Note: as of the `mc_backend/schema_finalize` branch, the Django project lives at the top-level `django/` directory, not `website/django`. Python dependencies are still declared in `website/requirements.txt`.

* **Install PostgreSQL 18** (matches production; see `CLAUDE.md`)
   * macOS: `brew install postgresql@18 && brew services start postgresql@18`. It's keg-only, so its binaries won't be on `PATH` by default — either symlink them or call them with the full path (`/opt/homebrew/opt/postgresql@18/bin/...`).
   * Linux: install from the PGDG apt repo.
* **Create the `postgres` role and `storage.pg` database** that `django/brca/site_settings.py` expects:
   ```
   createuser -s postgres
   psql -d postgres -c "ALTER USER postgres PASSWORD 'postgres';"
   createdb -O postgres storage.pg
   ```
* **Install the python dependencies**
   * `cd website`
   * `python3.13 -m venv <path-to-venv>` (a venv *outside* the repo, e.g. `~/.venv/brcaexchange`, works well — `website/.gitignore` also ignores an in-repo `env*` dir if you prefer that)
   * `source <path-to-venv>/bin/activate`
   * `pip install -r requirements.txt`
   * The project uses `psycopg` (v3), not `psycopg2`. If it can't find a working `libpq` at import time (common on Apple Silicon, where Homebrew's `postgresql@18` is keg-only and not on the default linker path), install the self-contained binary wheel instead: `pip install "psycopg[binary]==3.3.3"`. This avoids needing `DYLD_LIBRARY_PATH`/`libpq-dev` entirely — recommended for local dev.
* **Run the initial migration to populate the database**
   * `cd ../django` (top-level, not `website/django`)
   * `python manage.py migrate`
* **Start the server**
   * `python manage.py runserver` — serves on [http://localhost:8000/](http://localhost:8000/), `DEBUG=True` by default.

### Point the frontend at a local backend
`js/config.js` just re-exports `window.config` (`module.exports = window.config;`) — editing it directly no longer does anything. The actual default lives in `website/page.template`:
```js
window.config = {
    captcha_key: '...',
    backend_url: 'https://brcaexchange.org/backend',
    baseurl: '/'
}
```
In a real deployment this gets overridden by an Apache SSI include (`config.js`, injected per-environment — see `deployment/site_settings/`). `webpack-dev-server` doesn't process SSI, so for local development, temporarily edit `backend_url` above to `http://localhost:8000` — **don't commit this change**. (`CORS_ALLOWED_ORIGINS` in `django/brca/settings.py` already permits `http://localhost:8080`.)

### Lint

Use `npm run lint` to run the lint rules. We lint with eslint and babel-eslint.

## How to add data to your database
On `mc_backend/schema_finalize`, variant data during development is normally staged in a `pipeline` schema of `storage.pg` on a dev machine (e.g. `brcadev`), separate from the `public` schema the website's `default` DB connection actually reads. Table names are identical between the two schemas.

* Dump the `pipeline` schema from the source machine, e.g.:
  ```
  ssh -t brcadev 'sudo -u postgres pg_dump -d storage.pg -n pipeline -F c -f /tmp/pipeline_schema.dump'
  ssh brcadev 'sudo chown <you>:<you> /tmp/pipeline_schema.dump'
  scp brcadev:/tmp/pipeline_schema.dump .
  ```
* Restoring it directly won't populate `public` (and the schemas can't just be swapped — `public` also holds `auth`/`users`/`admin`/`sessions` tables that aren't in the pipeline dump). After `manage.py migrate` has created the `public` schema tables, extract the dump as data-only SQL, exclude `django_migrations` (it'll conflict with the bookkeeping `migrate` already wrote), rewrite the `pipeline.` schema prefix to `public.`, and load it with `session_replication_role = replica` set for the session (Postgres enforces FKs as internal triggers, and the plain-text dump doesn't preserve dependency order, so this avoids ordering errors — mirrors what `pg_restore` normally handles for you automatically):
  ```
  pg_restore --data-only -f /tmp/data.sql pipeline_schema.dump
  # manually strip the django_migrations COPY block, then:
  sed 's/pipeline\./public./g' /tmp/data.sql > /tmp/data_public.sql
  psql -d storage.pg -c "SET session_replication_role = replica;" -c "\i /tmp/data_public.sql" -c "SET session_replication_role = DEFAULT;"
  ```
* A few sequences in the pipeline dump (e.g. `enigma_domain_id_seq`, `genomic_coordinates_id_seq`) use older, unprefixed names that don't match the `public` schema's Django-generated sequence names (`data_enigma_domain_id_seq`, `variant_genomic_coordinates_id_seq`) — check `\ds public.*` and fix up the `setval(...)` calls if `psql` errors on them.
* In-progress pipeline data may not yet satisfy every `NOT NULL` constraint the current schema declares (e.g. `variant_genomic_coordinates.end_pos`). For local dev only, it's reasonable to relax such a constraint directly on your local DB (`ALTER TABLE ... ALTER COLUMN ... DROP NOT NULL;`) rather than editing the migration/model.

For a full non-pipeline dump (older workflow, may not apply to this branch): `sudo -u postgres pg_dump -d {REMOTEDBNAME} -F c -c -f /PATH/TO/full_db.dump`, then `pg_restore /PATH/TO/full_db.dump -c -v -1 -d storage.pg`.

## How to add additional releases
This process will add an additional release to the database and rebuild the words table to update autocomplete.

**Unverified on `mc_backend/schema_finalize`** — the schema no longer has a `words` table, so the autocomplete-rebuild step this describes is likely out of date for this branch.

 * Obtain a release archive (should be of the format release-MM-DD-YY.tar.gz).
 * From the project root directory, run `./deployment/deploy-data local PATH/TO/release-MM-DD-YY.tar.gz`.

## How to add/approve new users on the community page
* Production: go to https://brcaexchange.org/backend/admin/ and follow necessary steps (the `/backend` prefix is added by the production reverse proxy).
* Locally: [http://localhost:8000/admin/](http://localhost:8000/admin/) (create an account first with `python manage.py createsuperuser`).

### References
 * http://blog.keithcirkel.co.uk/how-to-use-npm-as-a-build-tool/
 * http://webpack.github.io/
 * http://www.youtube.com/watch?v=VkTCL6Nqm6Y
