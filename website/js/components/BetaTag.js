'use strict';

import React from 'react';
import {OverlayTrigger, Popover} from "react-bootstrap";

export default class BetaTag extends React.Component {
    render() {
        const betaBadge = (
            <span className="badge" style={{cursor: this.props.hoverText && 'help', verticalAlign: this.props.verticalAlign, marginLeft: this.props.margin, backgroundColor: '#3f6fff'}}>
            Beta Feature
            </span>
        );

        if (this.props.hoverText) {
            const popper = (
                <Popover title="Beta Feature">
                    <span dangerouslySetInnerHTML={{__html: this.props.hoverText}} />
                </Popover>
            );

            return (
                <OverlayTrigger placement='right' overlay={popper}>
                {betaBadge}
                </OverlayTrigger>
            );
        }

        return betaBadge;
    }
}

BetaTag.defaultProps = {
    hoverText: null,
    margin: '0.5em',
    verticalAlign: 'text-bottom'
};
