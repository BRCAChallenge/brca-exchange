'use strict';

var React = require('react');
var PureRenderMixin = require('./PureRenderMixin'); // deep-equals version of PRM
var {Button} = require('react-bootstrap');

var Community = React.createClass({
    mixins: [PureRenderMixin],
    render: function () {
        return (<Button href="/signup">Join our mailing list and this community space</Button>);
    }
});

module.exports = Community;