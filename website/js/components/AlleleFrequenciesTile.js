/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";
import {Panel, Table} from 'react-bootstrap';
import util from '../util';
const _ = require('underscore');
const KeyInline = require('./KeyInline');


export default class AlleleFrequenciesTile extends React.Component {
    constructor(props) {
        super(props);

        this.state = {
            // reportExpanded: this.defaultReportExpansions()
        };

        // this.reportToggled = this.reportToggled.bind(this);
        // this.setAllReportExpansion = this.setAllReportExpansion.bind(this);
    }

    getRowsAndDetermineIfEmpty(source, data, variant) {
        let rowsEmpty = 0;

        const rows = _.map(data, (rowDescriptor) => {
            let {prop, title} = rowDescriptor;
            let rowItem;

            if (variant[prop] !== null) {
                rowItem = util.getFormattedFieldByProp(prop, variant);
            }

            let isEmptyValue = util.isEmptyField(variant[prop]);

            if (isEmptyValue) {
                rowsEmpty += 1;
                rowItem = '-';
            }

            return (
                <tr key={prop} className={ (isEmptyValue && this.props.hideEmptyItems) ? "variantfield-empty" : "" }>
                    { rowDescriptor.tableKey !== false &&
                        <KeyInline tableKey={title} onClick={(event) => this.props.showHelp(event, title)}/>
                    }
                    <td colSpan={rowDescriptor.tableKey === false ? 2 : null} ><span className={"row-value" }>{rowItem}</span></td>
                </tr>
            );
        });
        const allEmpty = rowsEmpty >= data.length;

        return [rows, allEmpty];
    }

    generateHeader(field, fieldName) {

        return (
            <div className={`allele-frequency-header ${this.props.expanded ? 'expanded' : ''}`} onClick={this.onHandleToggle}>
                <div className="allele-frequency-cell allele-frequency-label">
                    {
                        this.props.expanded
                            ? <i className="fa fa-caret-down" aria-hidden="true" />
                            : <i className="fa fa-caret-right" aria-hidden="true" />
                    }
                    &nbsp;
                    <span>{fieldName}</span>
                </div>

                <div
                    className={`allele-frequency-cell allele-frequency-name ${!this.props.expanded ? 'collapsed' : ''}`}
                    style={{textAlign: 'left', flex: '1 1 auto'}}
                >
                </div>
            </div>
        );
    }

    generateCategoryTable(field, fieldName) {
        let renderedField;
        if (Array.isArray(renderedField)) {
            renderedField = field[0];
        } else {
            renderedField = field;
        }

        // let styles = this.getCollapsableClassSet();
        // const {submitter, cols, data} = this.props;

        // for each panel, construct key-value pairs as a row of the table
        // const submitterRows = cols.map(({prop, title, value, helpKey}) => {
        //     const isEmptyValue = util.isEmptyField(value);
        //     const rowItem = util.getFormattedFieldByProp(prop, data);

        //     return (
        //         <tr key={prop} className={ (isEmptyValue && this.props.hideEmptyItems) ? "variantfield-empty" : "" }>
        //             {
        //                 helpKey
        //                     ? <KeyInline tableKey={title} onClick={(event) => this.props.showHelp(event, helpKey)}/>
        //                     : <td className='help-target'><span style={{fontWeight: 'bold'}}>{title}</span></td>
        //             }
        //             <td><span className="row-value">{rowItem}</span></td>
        //         </tr>
        //     );
        // });

        return (
            <div>
                <div style={{marginBottom: 0, borderTop: 'solid 2px #ccc'}}>
                {
                    this.generateHeader(field, fieldName)
                }
                </div>

                <div ref='panel'>
                    <Table key={`frequency-name-${fieldName}`}>
                        <tbody>
                            {renderedField}
                        </tbody>
                    </Table>
                </div>
            </div>
        );
    }

    render() {
        const variant = this.props.variant;
        const data = this.props.alleleFrequencyData;

        const exacGraph = _.find(data, function(dd) {
                                return dd.source === "ExAC";
                            }).chart[0];
        const renderedExacGraph = exacGraph.replace(variant, exacGraph.prop);

        const exacData = _.find(data, function(dd) {
                                return dd.source === "ExAC";
                            }).data;
        const renderedExacData = this.getRowsAndDetermineIfEmpty("ExAC", exacData, variant);

        const thousandGenomesGraph = _.find(data, function(dd) {
                                return dd.source === "1000 Genomes";
                            }).chart[0];
        const renderedThousandGenomesGraph = thousandGenomesGraph.replace(variant, thousandGenomesGraph.prop);

        const thousandGenomesData = _.find(data, function(dd) {
                                return dd.source === "1000 Genomes";
                            }).data;
        const renderedThousandGenomesData = this.getRowsAndDetermineIfEmpty("1000 Genomes", thousandGenomesData, variant);

        const espData = _.find(data, function(dd) {
                                return dd.source === "ESP";
                            }).data;
        const renderedEspData = this.getRowsAndDetermineIfEmpty("ESP", espData, variant);

        const allEmpty = _.every([renderedExacData, renderedThousandGenomesData, renderedEspData], function(data) {
                                return data[1] === true;
                            }) && !variant.Variant_in_1000_Genomes && !variant.Variant_in_ExAC;

        // create the source panel itself now
        const groupTitle = `source-panel-${this.props.sourceName}`;
        const header = (
            <h3 style={{display: 'flex', flexDirection: 'row'}}>
                <a style={{flexGrow: 1}} href="#" onClick={(event) => this.props.onChangeGroupVisibility(groupTitle, event)}>
                    {this.props.groupTitle}
                </a>

                <a title='collapse all reports'
                    onClick={(event) => this.setAllReportExpansion(event, false)}
                    style={{cursor: 'pointer', marginRight: '10px'}}>
                    <i className="fa fa-angle-double-up" aria-hidden="true" />
                </a>

                <a title='expand all reports'
                    onClick={(event) => this.setAllReportExpansion(event, true)}
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
                    defaultExpanded={localStorage.getItem("collapse-group_" + groupTitle) !== "true"}>
                    {this.generateCategoryTable(renderedExacGraph, 'ExAC Graph')}
                    {this.generateCategoryTable(renderedExacData, 'ExAC Data')}
                    {this.generateCategoryTable(renderedThousandGenomesGraph, '1000 Genomes Graph')}
                    {this.generateCategoryTable(renderedThousandGenomesData, '1000 Genomes Data')}
                    {this.generateCategoryTable(renderedEspData, 'ESP Data')}
                </Panel>
            </div>
        );
    };
}
