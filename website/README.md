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
   * `pip install -qU -r requirements.txt`, if you encounter failures with psycopg2, try `pip install -U -r requirements.txt` and diagnose errors accordingly. You may need to `pip install psycopg2==2.7.3.2`.
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

## How to add data to your database
This process will add data to your database.

* Contact another developer to obtain a dump of the data from the beta or production site. The dump can be obtained by logging into the desired machine and running `sudo -u postgres pg_dump -d {REMOTEDBNAME} -F c -c -f /PATH/TO/full_db.dump` against whatever database is used as the source.

* Once you have the data on your machine, run `pg_restore /PATH/TO/full_db.dump -c -v -1 -d storage.pg` to load the data into your local database.

## How to add additional releases
This process will add an additional release to the database and rebuild the words table to update autocomplete.

 * Obtain a release archive (should be of the format release-MM-DD-YY.tar.gz).
 * From the project root directory, run `./deployment/deploy-data local PATH/TO/release-MM-DD-YY.tar.gz`.

## How to add/approve new users on the community page
* Go to http://brcaexchange.org/backend/admin/ and follow necessary steps.

### References
 * http://blog.keithcirkel.co.uk/how-to-use-npm-as-a-build-tool/
 * http://webpack.github.io/
 * http://www.youtube.com/watch?v=VkTCL6Nqm6Y
