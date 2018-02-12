# website for BRCA challenge
To contribute:
Fork the repository, make the changes, and submit pull request.

## Build the frontend
The build is based on npm and webpack.
* **Ensure that git and node are installed**
   * On OSX, install brew http://brew.sh/
       * `brew install git`
       * `brew install node`
   * On Linux:
        * `sudo apt-get install git-all`
        * `curl -sL https://deb.nodesource.com/setup_5.x | sudo -E bash -`
        * `sudo apt-get install -y nodejs`
* **Start the frontend**
   * `cd website`
   * `npm install`
      * Note that you may encounter problems with node-gyp installing contextify. Try reinstalling after switching to Node 6.9.1. NVM is a recommended version manager for Node.
   * `npm start`

## Build the server
The server runs on Django with postgres so install and set up those
* **Install postgres**
   * For Mac:
       * The easiest way is to install postgress.app from http://postgresapp.com/
       * Alternatively, `brew install postgresql`
   * Linux:
       * `sudo apt-get install postgresql postgresql-contrib`
* **Create a database called storage.pg**
   * `sudo -u postgres createdb storage.pg`
   * If this doesn't work, figure out how to create a db called `storage.pg` for user `postgres`.
* **Set the postgres role's password**
   * `sudo -u postgres psql postgres`
   *  at the prompt type `\password postgres` to set the password to `postgres`
   * If this doesn't work, figure out how to set user `postgres`'s password to `postgres`.
* **Install the python dependencies**
   * `cd website`
   * `pip install -qU -r requirements.txt`, if you encounter failures with psycopg2, try `pip install psycopg2==2.7.3.2`
* **Run the initial migration to populate the database**
   * `cd django`
   * `python manage.py migrate`
* **Start the server**
   * `python manage.py runserver`
* If you're developing locally point the frontend at the local server
   * Edit databaseUrl in `js/backend.js` to point to `localhost:8000`
* **Browse to [http://localhost:8080/](http://localhost:8080/)**

### Lint

Use `npm run lint` to run the lint rules. We lint with eslint and babel-eslint.


## How to change the data file
This process will delete all the previous data from the `variants` table and replace it with data from a new tsv file.

 * Replace the contents of `data/resources/aggregated.tsv` with the new data file
 * (Optional step if the schema has changed) Change the following files to correspond to schema changes (renamed / added / removed columns)
    *  `django/data/models.py` - this is the python model that corresponds to all the table columns
    *  `django/data/migrations/0001_initial.py` - specifies all the table columns, it matches the model above and the tsv headers
    *  `django/data/migrations/0002_search_index.py` - specifies which columns are used in full text search
    *  (if applicable) `django/data/migrations/0004_autocomplete_words.py` - specifies which columns are used in autocomplete suggestions
    *  `js/VariantTable.js`:
        * `columns` specifies which columns appear in the default mode and their names
        * `research_mode_columns` specifies which columns appear in research mode and their names
        * `subColumns` specifies which columns appear in the column selection filter and how they're grouped.
 * `cd django`
 * `python manage.py migrate --fake data zero && python manage.py migrate`

## How to add additional releases
This process will add an additional release to the database and rebuild the words table to update autocomplete.

 * Create a new directory at `brca-exchange/website/django/data/resources/releases/release-{MM-DD-YY}` that contains `built_with_change_types.tsv`, `removed.tsv`, and `version.json` from the desired pipeline output.
 * Add the archive of the release (e.g. `release-10-06-16.tar.gz`) to `brca-exchange/website/django/uploads/releases/` so the archive can be downloaded by users.
 * Run `./manage.py addrelease ./data/resources/releases/release-{MM-DD-YY}/built_with_change_types.tsv ./data/resources/releases/release-{MM-DD-YY}/reports.tsv ./data/resources/releases/release-{MM-DD-YY}/version.json ./data/resources/releases/release-{MM-DD-YY}/removed.tsv ./data/resources/releases/release-{MM-DD-YY}/diff.json` to add new data to the database, rebuild the words table, and generate a new release.

## How to add/approve new users on the community page
* Go to http://brcaexchange.org/backend/admin/ and follow necessary steps.

### References
 * http://blog.keithcirkel.co.uk/how-to-use-npm-as-a-build-tool/
 * http://webpack.github.io/
 * http://www.youtube.com/watch?v=VkTCL6Nqm6Y
