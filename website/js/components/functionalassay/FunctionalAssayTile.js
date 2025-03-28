'use strict';
import React from 'react';
import d3 from "d3";

const _ = require('underscore');

import CollapsibleTile from "../collapsible/CollapsibleTile";
import CollapsibleSection from "../collapsible/CollapsibleSection";
import FuncClassSubtile from "./FuncClassSubtile";
import {Table} from "react-bootstrap";
const util = require('../../util');
const slugify = require('../../slugify');
var KeyInline = require('../KeyInline');
const FunctionalAssayConstants = require("./FunctionalAssayConstants");


export const impacts = [
    { range: [-Infinity, -1.328], label: "Loss of Function", color: '#fc0d1b' },
    { range: [-1.328, -0.748], label: "Uncertain", color: '#fee49d' },
    { range: [-0.748, Infinity], label: "No Functional Impact", color: '#0f0' }
];

export default class FunctionalAssayTile extends React.Component {
    generateHeader(result) {
        return (
            <div className="func-assay-extras">
                <span className="func-assay-result" style={{float: 'right', paddingRight: '10px'}}>Report: {result}</span>
                <div style={{clear: 'both'}}></div>
            </div>
        );
    }

    render() {
        const results = this.props.results;
        const submitters = this.props.innerGroups[0].submitters[0];
        const funcScoreFindlay = parseFloat(this.props.variant.Functional_Enrichment_Findlay_ENIGMA_BRCA12_Functional_Assays);

        let allEmpty = !Object.values(results).some(function(val) {
            return !util.isEmptyField(val);
        });

        const impactScale = d3.scale.threshold()
            .domain(impacts[1].range)
            .range(impacts);

        let sections = _.map(FunctionalAssayConstants, (assay) => {
            let result = results[assay.name] ? results[assay.name] : "-";
            if (!util.isEmptyField(result) || !this.props.hideEmptyItems) {
                const rows = submitters[assay.name].map(({prop, title}) => {
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

                const additionalRowData = [
                    {key: 'Author', value: assay.author},
                    {key: 'Publication', value: assay.publication},
                    {key: 'Previous Publications', value: assay.previousPublications}
                ];

                const additionalRows = additionalRowData.map( field => {
                    let k = field.key;
                    let v = field.value;
                    const isEmptyValue = util.isEmptyField(v);
                    if (isEmptyValue) {
                        v = '-';
                    } else {
                        if (k === "Publication") {
                            // grab just the pmid for display
                            let splitText = v.split('/');
                            v = (<a href={v} target='_blank'>PMID:{splitText[splitText.length - 1]}</a>);
                        } else if (k === "Previous Publications") {
                            let pmids = v.split(',');
                            let length = pmids.length;
                            let links = pmids.map( (pmid, idx) => {
                                pmid = pmid.split(':')[1];
                                if (idx < length - 1) {
                                    return (<span><a href={`https://www.ncbi.nlm.nih.gov/pubmed/${pmid}`} target='_blank'>PMID:{pmid}</a>, </span>);
                                }
                                return (<span><a href={`https://www.ncbi.nlm.nih.gov/pubmed/${pmid}`} target='_blank'>PMID:{pmid}</a></span>);
                            });
                            v = links;
                        }
                    }
                    return (
                        <tr key={k} className={ (isEmptyValue && this.props.hideEmptyItems) ? "variantfield-empty" : "" }>
                            <KeyInline tableKey={k} noHelpLink={false}
                                tooltip={this.props.tooltips && this.props.tooltips[slugify(k)]}
                                onClick={(event) => this.props.showHelp(event, k)}
                            />
                            <td><span className="row-value">{v}</span></td>
                        </tr>
                    );
                });

                const resultDescriptionRow = [{key: 'Report Descriptions', value: assay.resultDescription}].map( field => {
                    let k = field.key;
                    let v = field.value;
                    const isEmptyValue = util.isEmptyField(v);
                    if (isEmptyValue) {
                        v = '-';
                    }
                    return (
                        <tr key={k} className={ (isEmptyValue && this.props.hideEmptyItems) ? "variantfield-empty" : "" }>
                            <KeyInline tableKey={k} noHelpLink={false}
                                tooltip={this.props.tooltips && this.props.tooltips[slugify(k)]}
                                onClick={(event) => this.props.showHelp(event, k)}
                            />
                            <td><span className="row-value">{v}</span></td>
                        </tr>
                    );
                });

                const allRows = additionalRows.concat(rows).concat(resultDescriptionRow);

                if (assay.name === "Findlay") {
                    return ( <CollapsibleSection
                                fieldName={assay.author}
                                extraHeaderItems={this.generateHeader(result)}
                                twoColumnExtraHeader={true}
                                defaultVisible={false}
                                assay={assay}
                                id={assay.name}
                            >
                                <div>
                                    <div className="subtile-container" style={{padding: 20, paddingTop: 10, paddingBottom: 10}}>
                                        This histogram summarizes functional classification results from saturation genome editing, provided by the Starita lab at the University of Washington (<a href="https://www.ncbi.nlm.nih.gov/pubmed/30209399">Findlay et al, Nature, 2018</a>). Please see <a href="https://sge.gs.washington.edu/BRCA1/">their site</a> for further information.  Please also note unpublished interim report for updated results for variants in Exon 5 from <a href="https://www.cangene-canvaruk.org/_files/ugd/ed948a_0399a952a1dc4767bed4364a04f6408b.pdf">Dace and Findlay 2023</a>
                                    </div>

                                    <div id="func-assay-container" style={{textAlign: 'center'}}>
                                        <FuncClassSubtile score={funcScoreFindlay} impactScale={impactScale} />
                                    </div>

                                    <div style={{padding: 20, paddingTop: 10, paddingBottom: 10, fontSize: 'smaller', fontStyle: 'italic', textAlign: 'center'}}>
                                        The median synonymous SNV scored 0.0 and the median nonsense SNV scored -2.12.
                                    </div>
                                    <div>
                                        <Table key={`functional-assay-name-${assay.name}`} >
                                            <tbody>
                                                {allRows}
                                            </tbody>
                                        </Table>
                                    </div>
                                </div>
                            </CollapsibleSection> );
                } else {
                    return ( <CollapsibleSection
                                fieldName={assay.author}
                                extraHeaderItems={this.generateHeader(result)}
                                twoColumnExtraHeader={true}
                                defaultVisible={false}
                                assay={assay}
                                id={assay.name}
                            >
                                <div>
                                    <Table key={`functional-assay-name-${assay.name}`} >
                                        <tbody>
                                            {allRows}
                                        </tbody>
                                    </Table>
                                </div>
                            </CollapsibleSection> );
                }
            }
        });

        let filteredSections = sections.filter(function(section) {
             return section !== undefined;
        });

        return (
            <CollapsibleTile allEmpty={allEmpty} {...this.props}>
                <div className="tile-disclaimer">
                    Assays were selected by the ENIGMA consortium as high quality assays that met internal standards permitting estimation of sensitivity and specificity, using approaches as described in <a href="https://pubmed.ncbi.nlm.nih.gov/31131967/" target='_blank'>Parsons 2019</a>.
                    <ul>
                        <li>Genetic variation can impact RNA transcripts and/or protein function or abundance. Assays labeled with a <b>P</b> measure the impact of genetic variation at the protein level. Assays labeled with a <b>R/P</b> measure impact at both the RNA and protein level.</li>
                        <li><b>Reports</b>, or the authors’ interpretation of the functional impact of each variant, are presented as they appear in the publications. A report of ‘Many Provided’ means that the authors performed more than one assay for each variant but provided no overall report of functional impact.</li>
                        <li>Additional assays may be added in the future.</li>
                    </ul>
                </div>
                <Table style={{paddingBottom: 0, marginBottom: 0}}>
                    <tbody>
                        <tr>
                            <td style={{verticalAlign: 'middle', width: '25%', fontWeight: 'bold', textAlign: 'center'}}>Disclaimer</td>
                            <td style={{verticalAlign: 'middle'}}>Functional assay data reflect laboratory models of disease, and should not be used as a substitute for clinical variant interpretation.  Further, please note that irrespective of the authors’ terms to reflect impact on function, the results presented <b>do not reflect a final variant classification</b>.</td>
                        </tr>
                    </tbody>
                </Table>

                {filteredSections}

            </CollapsibleTile>
        );
    }
}
