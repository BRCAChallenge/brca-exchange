/*eslint-env browser */
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

        // Popover id must be a valid HTML id (tableKey may include spaces/symbols)
        const popoverId = `tooltip_${String(tableKey).replace(/[^a-zA-Z0-9_-]/g, "_")}`;

        const popper = (
	    <Popover id={popoverId}>
                <Popover.Header as="h3">{tableKey}</Popover.Header>
                <Popover.Body>
                    <span dangerouslySetInnerHTML={{__html: tooltip}} />
                </Popover.Body>
            </Popover>
        );

        return (
            <td className='help-target'>
                <OverlayTrigger placement="bottom" trigger={['hover', 'focus']} overlay={popper}>
                    <span
                        className="help-target-inline"
                        role="button"
                        tabIndex={0}
                        onClick={onClick}
                        onKeyDown={(e) => {
                            if (e.key === 'Enter' || e.key === ' ') {
                                e.preventDefault();
                                onClick(e);
                            }
                        }}
                    >
                        {tableKey}
                    </span>
		</OverlayTrigger>
            </td>
        );
    }
}

export default KeyInline;
