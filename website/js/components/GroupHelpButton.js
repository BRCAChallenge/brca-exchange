/*eslint-env browser */
/*global require: false */
'use strict';

import React from 'react';

const GroupHelpButton = React.createClass({
    render() {
        const {onClick} = this.props;
        return (
            <span role='button' onClick={onClick} aria-label="Help" style={this.props.style}
                className='panel-help-btn glyphicon glyphicon-question-sign'
            />
        );
    }
});

export default GroupHelpButton;
