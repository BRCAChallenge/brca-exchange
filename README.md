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
 * `cd brca-website`
 * `npm install`
 * `npm start`
 * browse to [http://localhost:8080/webpack-dev-server/](http://localhost:8080/webpack-dev-server/)

### Lint

Use `npm run lint` to run the lint rules. We lint with eslint and babel-eslint.

### References
 * http://blog.keithcirkel.co.uk/how-to-use-npm-as-a-build-tool/
 * http://webpack.github.io/
 * http://www.youtube.com/watch?v=VkTCL6Nqm6Y
