/*eslint-env browser */
/*global require: false */
'use strict';

var React = require('react');
require('bootstrap/dist/css/bootstrap.css');
var Rx = require('rx');
var vcf = require('vcf.js');
require('rx.binding');
require('rx-dom');

var Row = require('react-bootstrap/lib/Row');
var Col = require('react-bootstrap/lib/Col');



var Application = React.createClass({

    render: function () {
        return (<div>hello world</div>);
    }
});

var main = document.getElementById('main');
React.render(<Application />, main);
