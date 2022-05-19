'use strict';
import React from 'react';


import * as _ from 'lodash';
import CollapsibleTile from "../collapsible/CollapsibleTile";
import CollapsibleSection from "../collapsible/CollapsibleSection";
import {Table} from "react-bootstrap";
const util = require('../../util');
const slugify = require('../../slugify');
import {geneMeta} from '../../SplicingData.js';
import {variantInfo, overlaps, pairwise} from '../../Splicing.js';
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

    getCIDomainOverlap(ciDomains, variantSpan) {
        const sortedDomains = _.sortBy(
            ciDomains.map(ciDomain => ({
                id: ciDomain.id,
                span: {
                    // the +1 makes the starting bound exclusive
                    name: ciDomain.span.name,
                    start: Math.min(ciDomain.span.start, ciDomain.span.end) + 1,
                    end: Math.max(ciDomain.span.start, ciDomain.span.end),
                }
            })),
            (a) => a.span.start
        );

        const overlappingDomains = sortedDomains
            .map((domain, idx) => ({idx: idx, domain: domain}))
            .filter(({domain}) => overlaps([variantSpan[0], variantSpan[1] - 1], [domain.span.start, domain.span.end]));


        return overlappingDomains.map(od => od.domain.span.name);
    }

    checkIfOverlapsExon(exons, variantSpan) {
        const sortedExons = _.sortBy(
            exons.map(exon => ({
                id: exon.id,
                span: {
                    // the +1 makes the starting bound exclusive
                    start: Math.min(exon.span.start, exon.span.end) + 1,
                    end: Math.max(exon.span.start, exon.span.end),
                }
            })),
            (a) => a.span.start
        );

        // sanity check: verify that the variant actually overlaps the gene at all
        const geneSpan = [sortedExons[0].span.start, sortedExons[sortedExons.length - 1].span.end];
        if (!overlaps(variantSpan, geneSpan)) {
            return null;
        }

        // we need to use the sorted set so that 'start' and 'end' are consistent regardless of the gene's direction
        let sorted = _.flatten(
            pairwise(sortedExons).map(([prevExon]) => {
                return {type: 'exon', ...prevExon};
            })
        );

        // and stick on the last element
        sorted.push({type: 'exon', ...(_.last(sortedExons))});

        // identify if variant overlaps exon
        // (variants are defined inclusively in the starting coordinate and exclusively in the ending one)
        const overlappingExons = sorted
            .map((exon, idx) => ({idx: idx, exon: exon}))
            .filter(({exon}) => overlaps([variantSpan[0], variantSpan[1] - 1], [exon.span.start, exon.span.end]));

        return (overlappingExons.length > 0) ? true : false;
    }

    getVarLoc(variant) {
        let info = variantInfo(variant);
        const variantStart = variant.Pos | 0;
        const variantEnd = variantStart + info.changed + info.deleted + info.inserted;
        const meta = geneMeta[variant['Gene_Symbol']];
        const exons = _.toPairs(meta.exons).map(([name, span]) => ({id: parseInt(name.substr(4)), span}));
        const ciDomains = _.toPairs(meta.CIDomains["ENIGMA Consortium"].domains).map(([name, span]) => ({id: name, span}));
        const variantSpan = [variantStart, variantEnd];

        let varLoc = this.getCIDomainOverlap(ciDomains, variantSpan);

        // only assign varLoc if variant overlaps CI domain AND is located in an exon
        if (varLoc.length > 0) {
            if (this.checkIfOverlapsExon(exons, variantSpan)) {
                return varLoc.join(', ');
            }
        }

        return 'N/A';
    }

    render() {
        const innerGroups = this.props.innerGroups;
        const variant = this.props.variant;
        let allEmpty = true;
        let sections = _.map(['varType', 'varLoc', 'BayesDel', 'SpliceAI'], (group) => {

            if (group === 'varType') {
                return ( <CollapsibleSection
                        fieldName={group}
                        computationalPrediction={true}
                        extraHeaderItems={this.generateHeader(variant.priors[group])}
                        twoColumnExtraHeader={true}
                        defaultVisible={false}
                        id={group}
                    />);
            } else if (group === 'varLoc') {
                let varLoc = this.getVarLoc(variant);
                return ( <CollapsibleSection
                        fieldName={group}
                        computationalPrediction={true}
                        extraHeaderItems={this.generateHeader(varLoc)}
                        twoColumnExtraHeader={true}
                        defaultVisible={false}
                        id={group}
                    />);
            }

            let groupData = innerGroups.filter(obj => {return obj.source === group;})[0].data;
            let result = variant[groupData.filter(obj => {return obj.title.includes('Result');})[0].prop];
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
