/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";
import {CollapsableMixin, Table} from "react-bootstrap";
import classNames from 'classnames';
import util from '../util';
import KeyInline from './KeyInline';

const VariantSubmitter = React.createClass({
    mixins: [CollapsableMixin],

    getCollapsableDOMNode: function() {
        return this.refs.panel.getDOMNode();
    },

    getCollapsableDimensionValue: function() {
        return this.refs.panel.getDOMNode().scrollHeight;
    },

    onHandleToggle: function(e) {
        e.preventDefault();
        this.setState({expanded: !this.state.expanded});
    },

    // FIXME: copied from index.js; find a way to use this without aliasing
    truncateData: function(field) {
        const fieldsToTruncate = ["Genomic_Coordinate_hg38", "Genomic_Coordinate_hg37", "Genomic_Coordinate_hg36"];
        return fieldsToTruncate.indexOf(field) > -1;
    },

    render: function() {
        let styles = this.getCollapsableClassSet();

        // ----

        const {submitter, cols} = this.props;

        // for each panel, construct key-value pairs as a row of the table
        const submitterRows = cols.map(({prop, title, value}) => {
            const rowItem = value;
            const isEmptyValue = util.isEmptyField(value);

            return (
                <tr key={prop} className={ (isEmptyValue && this.props.hideEmptyItems) ? "variantfield-empty" : "" }>
                    <KeyInline tableKey={title} onClick={(event) => this.showHelp(event, title)} />
                    <td colSpan={2} ><span className={ this.truncateData(prop) ? "row-value-truncated" : "row-value" }>{rowItem}</span></td>
                    <td>&nbsp;</td>
                </tr>
            );
        });

        // if we're not visible, these columns will be shown (if available)
        const dateLast = cols.find(x => x.prop === 'date-last-updated');

        return (
            <div>
                <Table key={`submitter-name-${submitter}`} style={{marginBottom: 0, borderTop: 'solid 2px #ccc'}}>
                    <thead>
                        <tr key="submitter-header" className={`submitter-header ${this.state.expanded ? 'expanded' : ''}`} onClick={this.onHandleToggle}>
                            <td className='help-target'>
                                {
                                    this.state.expanded
                                        ? <i className="fa fa-caret-down" aria-hidden="true" />
                                        : <i className="fa fa-caret-right" aria-hidden="true" />
                                }
                                &nbsp;
                                <span>Submitter:</span>
                            </td>
                            <td colSpan={2}><b>{submitter}</b></td>
                            <td style={{textAlign: 'right'}}>
                                { !this.state.expanded && dateLast && <span>(date: { dateLast.value })</span>}
                            </td>
                        </tr>
                    </thead>

                    <tbody>

                    </tbody>
                </Table>

                <div ref='panel' className={classNames(styles)}>
                    <Table key={`submitter-name-${submitter}`}>
                        <tbody>
                        { submitterRows }
                        </tbody>
                    </Table>
                </div>
            </div>
        );
    }
});

export default VariantSubmitter;
