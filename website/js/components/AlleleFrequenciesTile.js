/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";
import {Panel} from 'react-bootstrap';
import {AlleleFrequencyField} from "./AlleleFrequencyField";
import {AlleleFreqProvEvidField} from "./AlleleFrequencyProvisionalEvidenceField"; //ETK
import GroupHelpButton from './GroupHelpButton';

const _ = require('underscore');

// Includes the Allele Frequency fields plus the Provisional Evidence.
const fieldsOfInterest = {
    'Provisional ACMG Variant Evidence': true,
    'gnomAD V2.1 Exomes, Non-Cancer (Graphical)': true,
    'gnomAD V2.1 Exomes, Non-Cancer (Numerical)': false,
    'gnomAD V3.1 Genomes, Non-Cancer (Graphical)': true,
    'gnomAD V3.1 Genomes, Non-Cancer (Numerical)': false,
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
        console.log(this.props); // ETK
        const variant = this.props.variant;
        const data = this.props.alleleFrequencyData;

        const gnomadExomeGraph = [_.find(data, function(dd) {
                                return dd.source === "GnomAD Exomes";
                            }).chart[0], "gnomAD V2.1 Exomes, Non-Cancer (Graphical)"];

        const gnomadExomeData = [_.find(data, function(dd) {
                                return dd.source === "GnomAD Exomes";
                            }).data, "gnomAD V2.1 Exomes, Non-Cancer (Numerical)"];

        const gnomadGenomeGraph = [_.find(data, function(dd) {
                                return dd.source === "GnomADv3 Genomes";
                            }).chart[0], "gnomAD V3.1 Genomes, Non-Cancer (Graphical)"];

        const gnomadGenomeData = [_.find(data, function(dd) {
                                return dd.source === "GnomADv3 Genomes";
                            }).data, "gnomAD V3.1 Genomes, Non-Cancer (Numerical)"];

        const alleleFrequencyFields = [gnomadExomeGraph, gnomadExomeData, gnomadGenomeGraph,
                                       gnomadGenomeData];

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

        let fieldNameACMG = 'Provisional ACMG Variant Evidence';
        let expandedACMG = this.state[fieldNameACMG];
        let provisionalEvidence = {
            codePopfreq: variant.Provisional_Evidence_Code_Popfreq,
            descriptionPopfreq: variant.Provisional_Evidence_Description_Popfreq
        };
        const renderedProvEvidenceField = (
            <AlleleFreqProvEvidField
                field="AFPE Field"
                fieldName={fieldNameACMG}
                onFieldToggled={this.fieldToggled}
                expanded={expandedACMG}
                provisionalEvidence={provisionalEvidence}
            />
        );

        // TODO: figure out how to determine if everything is empty even though variant is in GnomAD
        const allEmpty = !variant.Variant_in_GnomAD;

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
                                    The gnomAD data sets used by BRCA Exchange are the “non-cancer” subsets
                                    of these sources. Data from <a href="https://tcga-data.nci.nih.gov/docs/publications/tcga/about.html">TCGA</a>
                                    &nbsp;and other cancer cohorts are excluded to ensure that the frequencies used
                                    to assess pathogenicity represent individuals not affected by cancer.
                                </div>
                            </div>
                            {renderedProvEvidenceField}
                            {renderedAlleleFrequencyFields}
                        </Panel.Body>
                    </Panel.Collapse>
                </Panel>
            </div>
        );
    }
};
