'use strict';
import React from 'react';

const _ = require('underscore');

import CollapsibleTile from "../collapsible/CollapsibleTile";
import CollapsibleSection from "../collapsible/CollapsibleSection";
import {Table} from "react-bootstrap";
const util = require('../../util');
const slugify = require('../../slugify');
const ComputationalPredictionConstants = require("./ComputationalPredictionConstants");
var KeyInline = require('../KeyInline');


export default class ComputationalPredictionTile extends React.Component {
    generateHeader(result) {
        return (
            <div className="func-assay-extras">
                <span className="func-assay-result" style={{float: 'right', paddingRight: '10px'}}>Result: {result}</span>
                <div style={{clear: 'both'}}></div>
            </div>
        );
    }

    render() {
        const innerGroups = this.props.innerGroups;
        const variant = this.props.variant;
        let allEmpty = true;
        let sections = _.map(['varType', 'varLoc', 'BayesDel', 'SpliceAI'], (group) => {

            if (group === 'varType' || group === 'varLoc') {
                return ( <CollapsibleSection
                        fieldName={group}
                        computationalPrediction={true}
                        extraHeaderItems={this.generateHeader(variant.priors[group])}
                        twoColumnExtraHeader={true}
                        defaultVisible={false}
                        id={group}
                    />);
            }

            let groupData = innerGroups.filter(obj => {return obj.source === group;})[0].data;
            let result = variant[groupData.filter(obj => {return obj.title === 'Result';})[0].prop];
            if (util.isEmptyField(result)) {
                result = '-';
            } else {
                allEmpty = false;
            };
            if (!util.isEmptyField(result) || !this.props.hideEmptyItems) {
                const rows = groupData.map(({prop, title}) => {
                    const isEmptyValue = util.isEmptyField(this.props.variant[prop]);
                    const rowItem = util.getFormattedFieldByProp(prop, this.props.variant);
                    return (
                        <tr key={prop} className={ (isEmptyValue && this.props.hideEmptyItems) ? "variantfield-empty" : "" }>
                            <KeyInline tableKey={title} noHelpLink={false}
                                tooltip={this.props.tooltips && this.props.tooltips[slugify(prop)]}
                                onClick={(event) => this.props.showHelp(event, prop)}
                            />
                            <td><span className="row-value">{rowItem}</span></td>
                        </tr>
                    );
                });

                let additionalRows = [];

                const additionalRowData = ComputationalPredictionConstants.filter(obj => {return obj.Method === group;});

                for (let [k, v] of Object.entries(additionalRowData[0])) {
                    const isEmptyValue = util.isEmptyField(v);
                    if (isEmptyValue) {
                        v = '-';
                    } else {
                        if (k === "Publication") {
                            additionalRows.push(
                                <tr key={k} className={ (isEmptyValue && this.props.hideEmptyItems) ? "variantfield-empty" : "" }>
                                    <KeyInline tableKey={k} noHelpLink={false}
                                        tooltip={this.props.tooltips && this.props.tooltips[slugify(k)]}
                                        onClick={(event) => this.props.showHelp(event, k)}
                                    />
                                    <td><span><a href={`https://www.ncbi.nlm.nih.gov/pubmed/${v}`} target='_blank'>PMID: {v}</a></span></td>
                                </tr>
                            );
                        }
                    }
                    additionalRows.push(
                        <tr key={k} className={ (isEmptyValue && this.props.hideEmptyItems) ? "variantfield-empty" : "" }>
                            <KeyInline tableKey={k} noHelpLink={false}
                                tooltip={this.props.tooltips && this.props.tooltips[slugify(k)]}
                                onClick={(event) => this.props.showHelp(event, k)}
                            />
                            <td><span className="row-value">{v}</span></td>
                        </tr>
                    );
                }

                const allRows = additionalRows.concat(rows);

                return ( <CollapsibleSection
                            fieldName={group}
                            computationalPrediction={true}
                            extraHeaderItems={this.generateHeader(result)}
                            twoColumnExtraHeader={true}
                            defaultVisible={false}
                            id={group}
                        >
                            <div>
                                <Table key={`computational-prediction-name-${group}`} >
                                    <tbody>
                                        {allRows}
                                    </tbody>
                                </Table>
                            </div>
                        </CollapsibleSection> );
            }
        });

        let filteredSections = sections.filter(function(section) {
             return section !== undefined;
        });

        return (
            <CollapsibleTile allEmpty={allEmpty} {...this.props}>
                <div className="tile-disclaimer">
                    This tile displays computational predictions and other annotations relevant for assigning ACMG/AMP bioinformatic codes.  Additional computational predictions are available via the <em>in Silico</em> Prior Predictions tile.
                </div>
                {filteredSections}
            </CollapsibleTile>
        );
    }
}
