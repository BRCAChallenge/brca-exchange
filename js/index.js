/*eslint-env browser */
/*global require: false */
'use strict';

var React = require('react');
require('bootstrap/dist/css/bootstrap.css');

var Application = React.createClass({
	render: function () {
		return (<div>Hello World</div>);
	}
});

var main = document.getElementById('main');
React.render(<Application />, main);
