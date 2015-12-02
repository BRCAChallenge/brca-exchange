# website for BRCA challenge
To contribute:
Fork the repository, make the changes, and submit pull request. 

## Build
The build is based on npm and webpack.
 * Ensure that git and node are installed
   * On OSX, install brew http://brew.sh/
   * `brew install git`
   * `brew install node`
 * `git clone https://github.com/BD2KGenomics/brca-website.git`
 * Ensure that the brca database is available to the build
   * Browse to http://brcaexchange.cloudapp.net/variants
   * Click "Download"
   * `mv ~/Downloads/variants.tsv .`
 * `cd brca-website`
 * `npm install`
 * `npm start`
 * browse to [http://localhost:8080/webpack-dev-server/](http://localhost:8080/webpack-dev-server/)

## Build with local server

 * Copy \*.table and \*.sqlite files into brca-website/backend/applications/data/databases
 * Edit url in js/backend.js to point to localhost:8000
 * Download and extract web2py
 * `python web2py.py -f ~/brca-website/backend`

For now, database files can be found on beta, in /var/www/backend/beta/web2py/applications/data/databases.

### Lint

Use `npm run lint` to run the lint rules. We lint with eslint and babel-eslint.

### References
 * http://blog.keithcirkel.co.uk/how-to-use-npm-as-a-build-tool/
 * http://webpack.github.io/
 * http://www.youtube.com/watch?v=VkTCL6Nqm6Y
