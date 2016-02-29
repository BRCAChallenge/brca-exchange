/*global module: false, require: false */
'use strict';

var React = require('react');

var RawHTML = React.createClass({
    render: function() {
        var {html, ...otherProps} = this.props;
        return (
            <div className='markdown' {...otherProps} dangerouslySetInnerHTML={{__html: html}} />
        );
    }
});

module.exports = RawHTML;