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
* **Checkout the repo**
   * `git clone https://github.com/BD2KGenomics/brca-website.git`
* **Start the frontend**
   * `cd brca-website`
   * `npm install`
   * `npm start`

## Build the server
The server runs on Django with postgres so install and set up those
* **Install postgres** 
   * For Mac:
       * The easiest way is to install postgress.app from http://postgresapp.com/
   * Linux: 
       * `sudo apt-get install postgresql postgresql-contrib`
* **Create the database** 
   * `sudo -u postgres createdb  storage.pg`
* **Set the postgres role's password** 
   * `sudo -u postgres psql postgres`
   *  at the prompt type `\password postgres` to set the password to `postgres`
* **Install the python dependencies**
   * `pip install -qU -r requirements.txt`
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

### References
 * http://blog.keithcirkel.co.uk/how-to-use-npm-as-a-build-tool/
 * http://webpack.github.io/
 * http://www.youtube.com/watch?v=VkTCL6Nqm6Y
