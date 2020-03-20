/*eslint-env browser */
/*global require: false */
'use strict';

import React from "react";
import {OverlayTrigger, Popover} from "react-bootstrap";

const KeyInline = React.createClass({
    render() {
        const { onClick, tableKey, tooltip, noHelpLink} = this.props;

        if (noHelpLink || !tooltip) {
            return <td className='help-target'><b>{tableKey}</b></td>;
        }

        const popper = (
            <Popover id={`tooltip_${tableKey}`} title={tableKey}>
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
});

export default KeyInline;
