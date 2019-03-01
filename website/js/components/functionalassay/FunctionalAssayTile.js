'use strict';
import React from 'react';
import d3 from "d3";

import CollapsibleTile from "../collapsible/CollapsibleTile";
import CollapsibleSection from "../collapsible/CollapsibleSection";
import FuncClassSubtile from "./FuncClassSubtile";
import {Table} from "react-bootstrap";


export const impacts = [
    { range: [-Infinity, -1.328], label: "Loss of Function", color: '#fc0d1b' },
    { range: [-1.328, -0.748], label: "Uncertain", color: '#fee49d' },
    { range: [-0.748, Infinity], label: "No Functional Impact", color: '#0f0' }
];

export default class FunctionalAssayTile extends React.Component {
    generateHeader(funcScore, label) {
        return (
            <div className="func-assay-extras">
                <span className="func-assay-score" style={{paddingLeft: '10px'}}>Score: {funcScore}</span>
                <span className="func-assay-impact" style={{float: 'right', paddingRight: '10px'}}>{label}</span>
                <div style={{clear: 'both'}}></div>
            </div>
        );
    }

    render() {
        const funcScore = parseFloat(this.props.score);

        const impactScale = d3.scale.threshold()
            .domain(impacts[1].range)
            .range(impacts);

        if (funcScore === undefined || isNaN(funcScore)) {
            return (
                <CollapsibleTile allEmpty={true} {...this.props}>
                    <div style={{padding: '10px'}}>
                    No functional assay data exists for this variant.
                    </div>
                </CollapsibleTile>
            );
        }

        return (
            <CollapsibleTile allEmpty={false} {...this.props}>
                <Table style={{paddingBottom: 0, marginBottom: 0}}>
                    <tbody>
                        <tr>
                            <td style={{verticalAlign: 'middle', width: '25%', fontWeight: 'bold', textAlign: 'center'}}>Disclaimer</td>
                            <td style={{verticalAlign: 'middle'}}>Please note that these experiments may include variants that have not been observed in clinical contexts. Furthermore, functional assay results reflect laboratory models of genetic variation and disease. Thus, these results should not be used as a proxy for clinical variant interpretations.</td>
                        </tr>
                    </tbody>
                </Table>
                <CollapsibleSection
                    fieldName="Findlay Function Score"
                    extraHeaderItems={this.generateHeader(funcScore.toFixed(3), impactScale(funcScore).label)}
                    twoColumnExtraHeader={true}
                    defaultVisible={true}
                >
                    <div className="subtile-container" style={{padding: 20, paddingTop: 10, paddingBottom: 10}}>
                        This histogram summarizes functional classification results from saturation genome editing, provided by the Starita lab at the University of Washington (<a href="https://www.ncbi.nlm.nih.gov/pubmed/30209399">Findlay et al, Nature, 2018</a>). Please see <a href="https://sge.gs.washington.edu/BRCA1/">their site</a> for further information.
                    </div>

                    <div id="func-assay-container" style={{textAlign: 'center'}}>
                        <FuncClassSubtile score={funcScore} impactScale={impactScale} />
                    </div>

                    <div style={{padding: 20, paddingTop: 10, paddingBottom: 10, fontSize: 'smaller', fontStyle: 'italic', textAlign: 'center'}}>
                        The median synonymous SNV scored 0.0 and the median nonsense SNV scored -2.12.
                    </div>

                    <ol style={{padding: 20, paddingTop: 10}}>
                        Publications that report on the accuracy of this assay:
                        <li style={{marginLeft: 30}}><a href="https://www.ncbi.nlm.nih.gov/pubmed/30209399">Findlay et al. 2018</a></li>
                    </ol>
                </CollapsibleSection>
            </CollapsibleTile>
        );
    }
}
