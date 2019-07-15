/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";
import {Panel} from 'react-bootstrap';
import {AlleleFrequencyField} from "./AlleleFrequencyField";
import GroupHelpButton from './GroupHelpButton';

const _ = require('underscore');


const fieldsOfInterest = {
    'GnomAD Genome (Graphical)': true,
    'GnomAD Genome (Numerical)': false,
    'GnomAD Exome (Graphical)': true,
    'GnomAD Exome (Numerical)': false,
    'ExAC (Graphical)': true,
    'ExAC (Numerical)': false,
    '1000 Genomes (Graphical)': false,
    '1000 Genomes (Numerical)': false,
    'ESP (Numerical)': false
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
            this.props.onFrequencyFieldToggled(this.collapser.getCollapsableDOMNode());
        });
    }

    fieldToggled(fieldName) {
        this.setState({
            [fieldName]: !this.state[fieldName]
        }, () => {
            // causes the parent to perform a (delayed) reflow
            this.props.onFrequencyFieldToggled(this.collapser.getCollapsableDOMNode());
        });
    }

    render() {
        const variant = this.props.variant;
        const data = this.props.alleleFrequencyData;

        const gnomadGenomesGraph = [_.find(data, function(dd) {
                                return dd.source === "GnomAD";
                            }).chart[0], "GnomAD Genome (Graphical)"];

        const gnomadGenomesData = [_.find(data, function(dd) {
                                return dd.source === "GnomAD";
                            }).data, "GnomAD Genome (Numerical)"];

        const gnomadExomesGraph = [_.find(data, function(dd) {
                                return dd.source === "GnomAD";
                            }).chart[0], "GnomAD Exome (Graphical)"];

        const gnomadExomesData = [_.find(data, function(dd) {
                                return dd.source === "GnomAD";
                            }).data, "GnomAD Exome (Numerical)"];

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

        const alleleFrequencyFields = [gnomadGenomesGraph, gnomadExomesGraph, gnomadGenomesData,
                                       gnomadExomesData, exacGraph, exacData, thousandGenomesGraph,
                                       thousandGenomesData, espData];

        const renderedAlleleFrequencyFields = alleleFrequencyFields.map((field, idx) => {
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
                    tooltips={this.props.tooltips}
                />
            );
        });

        // TODO: figure out how to determine if everything is empty even though variant is in 10KG or ExAC
        const allEmpty = !variant.Variant_in_1000_Genomes && !variant.Variant_in_ExAC && !variant.Variant_in_GnomAD;

        // create the source panel itself now
        const groupTitle = `source-panel-${this.props.sourceName}`;
        const header = (
            <h3>
                <a className="title" href="#" onClick={(event) => this.props.onChangeGroupVisibility(groupTitle, event, this.collapser.getCollapsableDOMNode())}>
                    {this.props.groupTitle}
                </a>

                <a title='collapse all fields'
                   className="toggle-subfields"
                    onClick={(event) => this.setAllFieldsExpansion(event, false)}
                    style={{cursor: 'pointer', marginRight: '10px'}}>
                    <i className="fa fa-angle-double-up" aria-hidden="true" />
                </a>

                <a title='expand all fields'
                   className="toggle-subfields"
                    onClick={(event) => this.setAllFieldsExpansion(event, true)}
                    style={{cursor: 'pointer'}}>
                    <i className="fa fa-angle-double-down" aria-hidden="true" />
                </a>

                {
                    this.props.helpSection &&
                    <GroupHelpButton group={this.props.helpSection}
                        onClick={(event) => {
                            this.props.showHelp(event, this.props.helpSection);
                            return true;
                        }}
                    />
                }
            </h3>
        );

        return (
            <div key={`group_collection-${groupTitle}`} className={ allEmpty && this.props.hideEmptyItems ? "group-empty variant-detail-group" : "variant-detail-group" }>
                <Panel
                    ref={(me) => { this.collapser = me; }}
                    header={header}
                    collapsable={true}
                    defaultExpanded={localStorage.getItem("collapse-group_" + groupTitle) !== "true"}
                    hideEmptyItems={this.props.hideEmptyItems}>
                    {renderedAlleleFrequencyFields}
                </Panel>
            </div>
        );
    }
};
