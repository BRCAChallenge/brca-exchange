/*eslint-env browser */
/*global require: false */
'use strict';

import React from "react";
import {OverlayTrigger, Popover} from "react-bootstrap";

const KeyInline = React.createClass({
    render() {
        const { onClick, tableKey, tooltip} = this.props;

        if (tooltip) {
            const popper = (
                <Popover title={tableKey}>
                    <span dangerouslySetInnerHTML={{__html: tooltip}} />
                </Popover>
            );

            return (
                <td className='help-target'>
                    <OverlayTrigger placement='bottom' overlay={popper}>
                        <span className="help-target-inline" onClick={onClick}>{tableKey}</span>
                    </OverlayTrigger>
                </td>
            );
        }
        else {
            return (
                <td className='help-target'>
                    <span className="help-target-inline" onClick={onClick}>{tableKey}</span>
                </td>
            );
        }
    }
});

export default KeyInline;
