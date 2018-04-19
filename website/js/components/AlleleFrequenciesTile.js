/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";
import {Panel} from 'react-bootstrap';
import {AlleleFrequencyField} from "./AlleleFrequencyField";
const _ = require('underscore');


const fieldsOfInterest = {
    'ExAC (Graphical)': true,
    'ExAC (Numerical)': true,
    '1000 Genomes (Graphical)': true,
    '1000 Genomes (Numerical)': true,
    'ESP (Numerical)': true
};

export default class AlleleFrequenciesTile extends React.Component {
    constructor(props) {
        super(props);

        this.state = fieldsOfInterest;

        this.fieldToggled = this.fieldToggled.bind(this);
        this.setAllFieldsExpansion = this.setAllFieldsExpansion.bind(this);
    }

    setAllFieldsExpansion(e, newExpansion) {
        e.stopPropagation();

        this.setState(function() {
            // updates all values in fieldsOfInterest to match newExpansion
            let updatedFieldsOfInterest = _.clone(fieldsOfInterest);
            Object.keys(updatedFieldsOfInterest).forEach( function(key) {
                updatedFieldsOfInterest[key] = newExpansion;
            });
            return updatedFieldsOfInterest;
        }, () => {
            // causes the parent to perform a (delayed) reflow
            this.props.onFrequencyFieldToggled();
        });
    }

    fieldToggled(fieldName) {
        this.setState({
            [fieldName]: !this.state[fieldName]
        }, () => {
            // causes the parent to perform a (delayed) reflow
            this.props.onFrequencyFieldToggled();
        });
    }

    render() {
        const variant = this.props.variant;
        const data = this.props.alleleFrequencyData;

        const exacGraph = [_.find(data, function(dd) {
                                return dd.source === "ExAC";
                            }).chart[0], "ExAC (Graphical)"];

        const exacData = [_.find(data, function(dd) {
                                return dd.source === "ExAC";
                            }).data, "ExAC (Numerical)"];

        const thousandGenomesGraph = [_.find(data, function(dd) {
                                return dd.source === "1000 Genomes";
                            }).chart[0], "1000 Genomes (Graphical)"];

        const thousandGenomesData = [_.find(data, function(dd) {
                                return dd.source === "1000 Genomes";
                            }).data, "1000 Genomes (Numerical)"];

        const espData = [_.find(data, function(dd) {
                                return dd.source === "ESP";
                            }).data, "ESP (Numerical)"];

        const alleleFrequencyFields = [exacGraph, exacData, thousandGenomesGraph, thousandGenomesData, espData].map((field, idx) => {
            let fieldValue = field[0];
            let fieldName = field[1];
            let expanded = this.state[fieldName];
            return (
                <AlleleFrequencyField
                    key={idx}
                    idx={idx}
                    variant={variant}
                    field={fieldValue}
                    fieldName={fieldName}
                    expanded={expanded}
                    onFieldToggled={this.fieldToggled}
                    hideEmptyItems={this.props.hideEmptyItems}
                    showHelp={this.props.showHelp}
                />
            );
        });

        // TODO: figure out how to determine if everything is empty
        const allEmpty = false;
        // const allEmpty = _.every([renderedExacData, renderedThousandGenomesData, renderedEspData], function(data) {
        //                         return data[1] === true;
        //                     }) && !variant.Variant_in_1000_Genomes && !variant.Variant_in_ExAC;

        // create the source panel itself now
        const groupTitle = `source-panel-${this.props.sourceName}`;
        const header = (
            <h3 style={{display: 'flex', flexDirection: 'row'}}>
                <a style={{flexGrow: 1}} href="#" onClick={(event) => this.props.onChangeGroupVisibility(groupTitle, event)}>
                    {this.props.groupTitle}
                </a>

                <a title='collapse all fields'
                    onClick={(event) => this.setAllFieldsExpansion(event, false)}
                    style={{cursor: 'pointer', marginRight: '10px'}}>
                    <i className="fa fa-angle-double-up" aria-hidden="true" />
                </a>

                <a title='expand all fields'
                    onClick={(event) => this.setAllFieldsExpansion(event, true)}
                    style={{cursor: 'pointer'}}>
                    <i className="fa fa-angle-double-down" aria-hidden="true" />
                </a>
            </h3>
        );

        return (
            <div key={`group_collection-${groupTitle}`} className={ allEmpty && this.state.hideEmptyItems ? "group-empty variant-detail-group" : "variant-detail-group" }>
                <Panel
                    header={header}
                    collapsable={true}
                    defaultExpanded={localStorage.getItem("collapse-group_" + groupTitle) !== "true"}
                    hideEmptyItems={this.state.hideEmptyItems}>
                    {alleleFrequencyFields}
                </Panel>
            </div>
        );
    }
};

//         const allEmpty = _.every([renderedExacData, renderedThousandGenomesData, renderedEspData], function(data) {
//                                 return data[1] === true;
//                             }) && !variant.Variant_in_1000_Genomes && !variant.Variant_in_ExAC;
