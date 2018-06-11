/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";
import {CollapsableMixin} from "react-bootstrap";
import classNames from 'classnames';
// import util from '../../util';
// import KeyInline from '../KeyInline';
// const _ = require('underscore');


const CollapsibleSection = React.createClass({
    mixins: [CollapsableMixin],

    getCollapsableDOMNode: function() {
        return this.refs.panel.getDOMNode();
    },

    getCollapsableDimensionValue: function() {
        return this.refs.panel.getDOMNode().scrollHeight;
    },

    handleToggle: function(e, fieldName) {
        e.preventDefault();

        // ask our parent to toggle us
        this.props.onFieldToggled(fieldName);
    },

    generateHeader: function(field, fieldName, extraHeaderItems) {
        return (
            <div className={`allele-frequency-header ${this.props.expanded ? 'expanded' : ''}`} onClick={(e) => this.handleToggle(e, fieldName)}>
                <div className="allele-frequency-cell allele-frequency-label">
                    {
                        this.props.expanded
                            ? <i className="fa fa-caret-down" aria-hidden="true" />
                            : <i className="fa fa-caret-right" aria-hidden="true" />
                    }
                    &nbsp;
                    <span>{fieldName}</span>
                </div>

                {
                    extraHeaderItems &&
                    <div className="submitter-cell optional" style={{textAlign: 'left', flex: '0 1 auto'}}>
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
        const {field, fieldName, hideEmptyItems, extraHeaderItems} = this.props;
        let allEmpty = false;
        let styles = this.getCollapsableClassSet();

        return (
            <div className={ allEmpty && hideEmptyItems ? "group-empty" : "" }>
                <div style={{marginBottom: 0, borderTop: 'solid 2px #ccc'}}>
                { this.generateHeader(field, fieldName, extraHeaderItems) }
                </div>

                <div ref='panel' className={allEmpty ? "group-empty" : classNames(styles)}>
                {this.props.children}
                </div>
            </div>
        );
    }
});

CollapsibleSection.defaultProps = {
    defaultVisible: false
};

module.exports = CollapsibleSection;
export default CollapsibleSection;
