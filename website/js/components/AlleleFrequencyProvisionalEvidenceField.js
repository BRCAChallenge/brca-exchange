/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";
import ReactDOM from 'react-dom';
import {Collapse, Table} from "react-bootstrap";
import KeyInline from './KeyInline';


const AlleleFreqProvEvidField = React.createClass({

    getInitialState: function () {
        // identifies which subpopulation groups are expanded
        return {};
    },

    getCollapsableDOMNode: function() {
        return ReactDOM.findDOMNode(this.refs.panel);
    },

    getCollapsableDimensionValue: function() {
        return ReactDOM.findDOMNode(this.refs.panel).scrollHeight;
    },

    handleToggle: function(e, fieldName) {
        e.preventDefault();

        // ask our parent to toggle us
        this.props.onFieldToggled(fieldName);
    },

    fieldToggled: function(field) {
        // handles toggling of subpopulation groups in gnomad numerical tables
        this.setState({
            [field]: !this.state[field]
        });
    },

    generateHeader: function(field, fieldName) {
        console.log("ETK rendered acmg header for ." + field + "." + fieldName + ". expanded:" + this.props.expanded + ".");
        return (
            <div className={`allele-frequency-header ${this.props.expanded ? 'expanded' : ''}`} onClick={
                (e) => this.handleToggle(e, fieldName)}>
                <div className="allele-frequency-cell allele-frequency-label">
                    {
                        this.props.expanded
                            ? <i className="fa fa-caret-down" aria-hidden="true" />
                            : <i className="fa fa-caret-right" aria-hidden="true" />
                    }
                    &nbsp;
                    <span>{fieldName}</span>
                </div>
            </div>
        );
    },
    render: function() {
        const {field, fieldName} = this.props;
        return (
            <div className="allele-frequency-provisional-evidence">
                <div style={{marginBottom: 0, borderTop: 'solid 2px #ccc'}}>
                {
                    this.generateHeader(field, fieldName)
                }
                </div>

                <Collapse className="allele-frequency-provisional-evidence-collapse"
                    in={this.props.expanded}
                    onEntered={this.props.relayoutGrid}
                    onExited={this.props.relayoutGrid}
                >
                    <div>
                        <div className="tile-disclaimer">
                            <div>
                                The following describes the population frequency evidence under the Variant Curation
Expert Panel (VCEP) rules of the ENIGMA BRCA1/2 VCEP, Version 1.1.0, dated 2023-11-22.  <b>This provisional evidence
assessment has not yet been reviewed by the VCEP biocuration team.</b>
                            </div>
                        </div>
                        <Table>
                            <tbody>
                                <tr>
                                    {
                                        <KeyInline tableKey="Provisional Code" />
                                    }
                                    <td>{this.props.provisionalEvidence.codePopfreq}</td>
                                </tr>
                                <tr>
                                    {
                                        <KeyInline tableKey="Description" />
                                    }
                                    <td>{this.props.provisionalEvidence.descriptionPopfreq}</td>
                                </tr>


                            </tbody>
                        </Table>
                    </div>
                </Collapse>
            </div>
        );
    }
});

module.exports = {
    AlleleFreqProvEvidField,
};
