/*eslint-env browser */
/*global require: false */
'use strict';

import React from "react";

const KeyInline = React.createClass({
    render() {
        const {onClick, tableKey} = this.props;
        return (
            <td className='help-target'>
                <span className="help-target-inline" onClick={onClick}>{tableKey}</span>
            </td>
        );
    }
});

export default KeyInline;
