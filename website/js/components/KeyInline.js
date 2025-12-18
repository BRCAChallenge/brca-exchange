/*eslint-env browser */
/*global require: false */
'use strict';

import React from "react";
import {OverlayTrigger, Popover} from "react-bootstrap";

class KeyInline extends React.PureComponent {
    getCaret = () => {
        return (
            this.props.expanded
                ? <i className="fa fa-caret-down gnomad-header-row" aria-hidden="true" />
                : <i className="fa fa-caret-right gnomad-header-row" aria-hidden="true" />
        );
    };

    render() {
        const { onClick, tableKey, tooltip, noHelpLink, headerGroup } = this.props;

        if (noHelpLink || !tooltip) {
            return <td className='help-target'>{headerGroup ? this.getCaret() : ''}<b>{tableKey}</b></td>;
        }

        const popper = (
	    <Popover id={`tooltip_${tableKey}`}>
                <Popover.Header as="h3">{tableKey}</Popover.Header>
                <Popover.Body>
                    <span dangerouslySetInnerHTML={{__html: tooltip}} />
                </Popover.Body>
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
}

export default KeyInline;
