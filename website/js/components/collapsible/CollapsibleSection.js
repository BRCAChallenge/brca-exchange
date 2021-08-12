/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";
import {Collapse} from "react-bootstrap";

const CollapsibleSection = React.createClass({
    handleToggle: function(e, id) {
        e.preventDefault();

        // ask our parent to toggle us
        this.props.onFieldToggled(id);
    },

    generateHeader: function(id, fieldName, extraHeaderItems, twoColumnExtraHeader, assay) {
        if (assay) {
            let className = `${assay.loc.toLowerCase()}-header`;
            var assayLoc = (
                <span className="functional-assay-header"><span className={className}>{assay.loc}</span></span>
            );
        }

        return (
            <div className={`allele-frequency-header ${this.props.expanded ? 'expanded' : ''}`} onClick={(e) => this.handleToggle(e, id)}>
                <div className="allele-frequency-cell allele-frequency-label">
                    {assay ? assayLoc : ''}
                    {
                        this.props.expanded
                            ? <i className="fa fa-caret-down" aria-hidden="true" />
                            : <i className="fa fa-caret-right" aria-hidden="true" />
                    }
                    &nbsp;
                    <span>{fieldName}</span>
                </div>
                {
                    extraHeaderItems && !twoColumnExtraHeader &&
                    <div className="submitter-cell optional" style={{textAlign: 'left', flex: '0 1 auto'}}>
                        {
                            // remaining header elements depend on the source
                            extraHeaderItems
                        }
                    </div>
                }
                {
                    extraHeaderItems && twoColumnExtraHeader &&
                    <div className="submitter-cell optional" style={{textAlign: 'left', width: '100%'}}>
                        {
                            // remaining header elements depend on the source
                            extraHeaderItems
                        }
                    </div>
                }
            </div>
        );
    },

    render: function() {
        const {id, fieldName, hideEmptyItems, extraHeaderItems, twoColumnExtraHeader, assay} = this.props;

        let allEmpty = false;

        return (
            <div className={ allEmpty && hideEmptyItems ? "group-empty" : "" }>
                <div style={{marginBottom: 0, borderTop: 'solid 2px #ccc'}}>
                { this.generateHeader(id, fieldName, extraHeaderItems, twoColumnExtraHeader, assay) }
                </div>

                <Collapse className={allEmpty ? "group-empty" : ""}
                    in={this.props.expanded}
                    onEntered={this.props.relayoutGrid}
                    onExited={this.props.relayoutGrid}
                >
                    <div>
                    {this.props.children}
                    </div>
                </Collapse>
            </div>
        );
    }
});

CollapsibleSection.defaultProps = {
    defaultVisible: false
};

module.exports = CollapsibleSection;
export default CollapsibleSection;
