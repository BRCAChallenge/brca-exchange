/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";
import {Panel} from 'react-bootstrap';
import {AlleleFrequencyField} from "./AlleleFrequencyField";
import GroupHelpButton from './GroupHelpButton';

const _ = require('underscore');


const fieldsOfInterest = {
    'gnomAD Exomes (Graphical)': true,
    'gnomAD Exomes (Numerical)': false,
    'gnomAD Genomes (Graphical)': true,
    'gnomAD Genomes (Numerical)': false,
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
        });
    }

    fieldToggled(fieldName) {
        this.setState({
            [fieldName]: !this.state[fieldName]
        });
    }

    render() {
        const variant = this.props.variant;
        const data = this.props.alleleFrequencyData;

        const gnomadExomeGraph = [_.find(data, function(dd) {
                                return dd.source === "GnomAD Exomes";
                            }).chart[0], "gnomAD Exomes (Graphical)"];

        const gnomadExomeData = [_.find(data, function(dd) {
                                return dd.source === "GnomAD Exomes";
                            }).data, "gnomAD Exomes (Numerical)"];

        const gnomadGenomeGraph = [_.find(data, function(dd) {
                                return dd.source === "GnomAD Genomes";
                            }).chart[0], "gnomAD Genomes (Graphical)"];

        const gnomadGenomeData = [_.find(data, function(dd) {
                                return dd.source === "GnomAD Genomes";
                            }).data, "gnomAD Genomes (Numerical)"];

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


        const alleleFrequencyFields = [gnomadExomeGraph, gnomadExomeData, gnomadGenomeGraph,
                                       gnomadGenomeData, exacGraph, exacData, thousandGenomesGraph,
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
                    relayoutGrid={this.props.relayoutGrid}
                    showHelp={this.props.showHelp}
                    tooltips={this.props.tooltips}
                />
            );
        });

        // TODO: figure out how to determine if everything is empty even though variant is in 10KG / ExAC / GnomAD
        const allEmpty = !variant.Variant_in_1000_Genomes && !variant.Variant_in_ExAC && !variant.Variant_in_GnomAD;

        // create the source panel itself now
        const groupTitle = `source-panel-${this.props.sourceName}`;

        // return <div />;

        return (
            <div key={`group_collection-${groupTitle}`} className={ allEmpty && this.props.hideEmptyItems ? "group-empty variant-detail-group" : "variant-detail-group" }>
                <Panel
                    collapsible={true}
                    defaultExpanded={localStorage.getItem("collapse-group_" + groupTitle) !== "true"}
                    hideEmptyItems={this.props.hideEmptyItems}
                >
                    <Panel.Heading>
                        <Panel.Title componentClass="h3">
                            <Panel.Toggle componentClass="a" className="title"
                                onClick={(event) => this.props.onChangeGroupVisibility(groupTitle, event)}
                            >
                                {this.props.groupTitle}
                            </Panel.Toggle>

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
                        </Panel.Title>
                    </Panel.Heading>
                    <Panel.Collapse
                        onEntered={this.props.relayoutGrid}
                        onExited={this.props.relayoutGrid}
                    >
                        <Panel.Body>
                            <div className="tile-disclaimer">
                                <div>
                                    The gnomAD and ExAC data sets used by BRCA Exchange are the “non-cancer” subsets
                                    of these sources. Data from <a href="https://tcga-data.nci.nih.gov/docs/publications/tcga/about.html">TCGA</a>
                                    &nbsp;and other cancer cohorts are excluded to ensure that the frequencies used
                                    to assess pathogenicity represent individuals not affected by cancer.
                                </div>
                            </div>
                            {renderedAlleleFrequencyFields}
                        </Panel.Body>
                    </Panel.Collapse>
                </Panel>
            </div>
        );
    }
};
