/*global document: false, require: false, __webpack_public_path__: true */
'use strict';

// A minimalist router, using webpack to build chunks for each page.
// To add new pages, extend the condition with additional url patterns
// and require.ensure() blocks.

var config = require('config');

__webpack_public_path__ = config.baseurl; // eslint-disable-line camelcase

// XXX There is an uglify2 bug which will drop the require() calls if we don't access a
// property on the return value.
// https://github.com/mishoo/UglifyJS2/commit/276b9a31cda2a2ef93e7af4e966baae91a434361

var path = document.location.pathname;

if (path.match(/^\/$/)) {
	require.ensure(['index'], function () {
        // XXX see above
		require(['index']).foo; // eslint-disable-line no-unused-expressions
	});
} else {
    document.write('Page not found');
}
