/*global module: false, require: false */

import React from 'react';

var RawHTML = React.createClass({
    render: function() {
        var {html, ...otherProps} = this.props;
        return (
            <div className='markdown' {...otherProps} dangerouslySetInnerHTML={{__html: html}} />
        );
    }
});

export default RawHTML;
