/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";
import ReactDOM from 'react-dom';
import {Collapse, Table} from "react-bootstrap";
import classNames from "classnames";
import util from '../util';
import KeyInline from './KeyInline';
import slugify from "../slugify";
const _ = require('underscore');


const AlleleFrequencyField = React.createClass({

    getInitialState: function () {
        // identifies which subpopulation groups are expanded
        return {
            'Allele_frequency_exome_NFE_GnomAD': false,
            'Allele_frequency_exome_EAS_GnomAD': false,
            'Allele_frequency_genome_NFE_GnomAD': false,
            'Allele_frequency_genome_EAS_GnomAD': false
        };
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

    generateHeader: function(field, fieldName, flag) {
        let fnLower = fieldName.toLowerCase();
        let isGenome = false;
        let isGnomad = false;
        if (fnLower.includes('gnomad')) {
            isGnomad = true;
            if (fnLower.includes('genome')) {
                isGenome = true;
            }
        }

        return (
            <div className={`allele-frequency-header ${this.props.expanded ? 'expanded' : ''}`} onClick={(e) => this.handleToggle(e, fieldName)}>
                <div className="allele-frequency-cell allele-frequency-label">
                    {
                        isGnomad
                            ? isGenome
                                ? <span className="allele-frequency-gnomad-header"><span className="genome-header">G</span><span className="glyphicon glyphicon-flag" style={{display: flag && !util.isEmptyField(flag) ? '' : 'none'}}></span></span>
                                : <span className="allele-frequency-gnomad-header"><span className="exome-header">E</span><span className="glyphicon glyphicon-flag" style={{display: flag && !util.isEmptyField(flag) ? '' : 'none'}}></span></span>
                            : ''
                    }
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

    getPopMax: function(fieldName, variant) {
        /*eslint-disable camelcase*/
        let popmax = parseFloat(variant.faf95_popmax_genome_GnomADv3) ? parseFloat(variant.faf95_popmax_genome_GnomADv3).toPrecision(4) : '-';
        if (fieldName.includes("Genomes, Non-Cancer (Graphical)")) {
            return popmax + ' (' + variant.faf95_popmax_population_genome_GnomADv3 + ')';
        } else {
            return popmax + ' (' + variant.faf95_popmax_population_exome_GnomAD + ')';
        }
        /*eslint-enable camelcase*/
    },

    cleanRowTitle: function (title) {
        title = title.replace(" (gnomAD V2.1 Exomes)", "");
        title = title.replace(" (gnomAD V3.1 Genomes)", "");
        return title;
    },

    getRowsAndDetermineIfEmpty: function(source, data, variant) {
        let rowsEmpty = 0;

        const rows = _.map(data, (rowDescriptor) => {
            let {prop, title, noHelpLink} = rowDescriptor;
            let rowItem;

            title = this.cleanRowTitle(title);

            if (variant[prop] !== null) {
                rowItem = util.getFormattedFieldByProp(prop, variant);
            }

            let isEmptyValue = util.isEmptyField(variant[prop]);

            if (isEmptyValue) {
                rowsEmpty += 1;
                rowItem = '-';
            }

            let rowClasses = classNames({
                'variantfield-empty': (isEmptyValue && this.props.hideEmptyItems),
            });

            return (
                <tr key={prop} className={ rowClasses }>
                    { rowDescriptor.tableKey !== false &&
                        (
                            <KeyInline tableKey={title} noHelpLink={noHelpLink}
                                tooltip={this.props.tooltips && this.props.tooltips[slugify(prop)]}
                                onClick={(event) => this.props.showHelp(event, prop)}
                            />
                        )
                    }
                    <td colSpan={rowDescriptor.tableKey === false ? 2 : null} ><span className={"row-value" }>{rowItem}</span></td>
                </tr>
            );
        });
        const allEmpty = rowsEmpty >= data.length;
        return [rows, allEmpty];
    },

    render: function() {
        const {field, fieldName, variant, hideEmptyItems} = this.props;
        let renderedRows, flag, gnomadLink, chr, transcript, refSeqTranscript, lrgLink, variantId, dataset;
        let allEmpty = false;
        let isChart = false;
        let isGnomad = false;
        let popmax = this.getPopMax(fieldName, variant);


        if (fieldName.toLowerCase().includes('gnomad')) {
            flag = variant.Flags_GnomAD;
            chr = variant.Chr;
            isGnomad = true;
            if (fieldName.toLowerCase().includes('genome')) {
                variantId = variant.Variant_id_GnomADv3;
                dataset = "gnomad_r3_non_cancer";
            } else {
                variantId = variant.Variant_id_GnomAD;
                dataset = "gnomad_r2_1_non_cancer";
            }
            gnomadLink = "https://gnomad.broadinstitute.org/variant/" + variantId + "?dataset=" + dataset;
            if (chr === "17") {
                transcript = "ENST00000357654";
                refSeqTranscript = "NM_007294.3";
                lrgLink = "http://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_292.xml";
            } else {
                transcript = "ENST00000544455";
                refSeqTranscript = "NM_000059.3";
                lrgLink = "http://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_293.xml";
            }
        }

        if (fieldName === "gnomAD V3.1 Genomes, Non-Cancer (Graphical)") {
            renderedRows = field.replace(variant, field.prop);
            if (!variant.Variant_in_GnomAD || util.isEmptyField(variant.Allele_frequency_genome_GnomADv3)) {
                allEmpty = true;
            }
            isChart = true;
        } else if (fieldName === "gnomAD V3.1 Genomes, Non-Cancer (Numerical)") {
            renderedRows = this.getRowsAndDetermineIfEmpty("GnomAD", field, variant, flag);
        } else if (fieldName === "gnomAD V2.1 Exomes, Non-Cancer (Graphical)") {
            renderedRows = field.replace(variant, field.prop);
            if (!variant.Variant_in_GnomAD || util.isEmptyField(variant['Allele_frequency_exome_GnomAD'])) {
                allEmpty = true;
            }
            isChart = true;
        } else if (fieldName === "gnomAD V2.1 Exomes, Non-Cancer (Numerical)") {
            renderedRows = this.getRowsAndDetermineIfEmpty("GnomAD", field, variant, flag);
        } else if (fieldName === "ExAC (Graphical)") {
            renderedRows = field.replace(variant, field.prop);
            if (!variant.Variant_in_ExAC) {
                allEmpty = true;
            }
            isChart = true;
        } else if (fieldName === "ExAC (Numerical)") {
            renderedRows = this.getRowsAndDetermineIfEmpty("ExAC", field, variant);
        } else if (fieldName === "1000 Genomes (Graphical)") {
            renderedRows = field.replace(variant, field.prop);
            if (!variant.Variant_in_1000_Genomes) {
                allEmpty = true;
            }
            isChart = true;
        } else if (fieldName === "1000 Genomes (Numerical)") {
            renderedRows = this.getRowsAndDetermineIfEmpty("1000 Genomes", field, variant);
        } else if (fieldName === "ESP (Numerical)") {
            renderedRows = this.getRowsAndDetermineIfEmpty("ESP", field, variant);
        }

        if (Array.isArray(renderedRows)) {
            allEmpty = renderedRows[1];
            renderedRows = renderedRows[0];
        }

        return (
            <div className={ allEmpty && (isChart || hideEmptyItems) ? "group-empty" : "" }>
                <div style={{marginBottom: 0, borderTop: 'solid 2px #ccc'}}>
                {
                    this.generateHeader(field, fieldName, flag)
                }
                </div>

                <Collapse className={allEmpty && isChart ? "group-empty" : ""}
                    in={this.props.expanded}
                    onEntered={this.props.relayoutGrid}
                    onExited={this.props.relayoutGrid}
                >
                    <div>
                        {flag && flag !== '-'
                            ? <div className="glyphicon glyphicon-flag gnomad-flag"><span style={{color: 'black', marginLeft: '6px'}}>{flag}</span></div>
                            : ''
                        }
                        {flag && flag !== '-'
                            ? <hr style={{marginTop: "2px"}}/>
                            : ''
                        }
                        <div className="tile-disclaimer" style={{display: isGnomad && !util.isEmptyField(flag) ? '' : 'none'}}>
                            <div>
                                You are viewing flags for ENSEMBL transcript {transcript}. This is not the canonical
                                transcript shown by default on gnomAD, but corresponds to RefSeq transcript {refSeqTranscript}
                                &nbsp;(per <a href={lrgLink}>LRG</a>).
                            </div>
                        </div>
                        <div className="tile-disclaimer" style={{display: isGnomad ? '' : 'none'}}>
                            <div>
                                Additional data for this variant, including detailed
                                populations, quality scores, and flags relative to other transcripts,
                                <a href={gnomadLink} target="_blank">&nbsp;are available at gnomAD</a>.
                            </div>
                        </div>
                        <Table key={`allele-frequency-name-${fieldName}`} >
                            <tbody>
                            { renderedRows }
                            </tbody>
                        </Table>
                        <div className="tile-disclaimer" style={{display: isGnomad && isChart && !util.isEmptyField(popmax) ? '' : 'none'}}>
                            <div>
                                Popmax Filtering AF (95% confidence): { popmax }
                            </div>
                            <div>
                                If the filter allele frequency of a variant is above the maximum credible population AF
                                for a condition of interest, then this variant should be filtered (ie not considered a
                                candidate causative variant).
                            </div>
                        </div>
                    </div>
                </Collapse>
            </div>
        );
    }
});

module.exports = {
    AlleleFrequencyField,
};
