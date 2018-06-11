/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";

import CollapsibleTile from "../collapsible/CollapsibleTile";
import CollapsibleSection from "../collapsible/CollapsibleSection";

import SplicingLevelSubtile from "./SplicingLevelSubtile";
import ProteinLevelSubtile from "./ProteinLevelSubtile";
import InSilicoPredSubtile from "./InSilicoPredSubtile";

const dummyData = {
    pathosProb: {
        value: 0.5,
        context: [
            { title: 'Variant Location', prop: 'var_loc', value: 'Exon' },
            { title: 'Variant Type', prop: 'var_type', value: 'Missense' }
        ]
    },
    proteinLevelEst: {
        value: 0.02
    },
    splicingLevelEst: {
        value: 0.64,
        comparisons: {
            wildType: {
                sequence: 'TTGgtaaaa',
                rawScore: 3.23,
                zScore: -3.293
            },
            variant: {
                sequence: 'TTGgaaaaa',
                rawScore: -4.95,
                zScore: -5.6
            }
        },
        pathosProb: {
            outside: 0.02,
            improved: 0.04,
            low: 0.04,
            moderate: 0.34,
            high: 0.97
        }
    }
};

export default class SilicoPredTile extends React.Component {
    constructor(props) {
        super(props);
    }

    render() {
        return (
            <CollapsibleTile {...this.props}>
                <CollapsibleSection
                    fieldName="In Silico Probability of Pathogenicity"
                    extraHeaderItems={dummyData.pathosProb.value}
                >
                    <InSilicoPredSubtile data={dummyData.pathosProb} />
                </CollapsibleSection>

                <CollapsibleSection
                    fieldName="Protein-level Estimation"
                    extraHeaderItems={dummyData.proteinLevelEst.value}
                    defaultVisible={false}
                >
                    <ProteinLevelSubtile probability={dummyData.proteinLevelEst.value} />
                </CollapsibleSection>

                <CollapsibleSection
                    fieldName="Splicing-level Estimation"
                    extraHeaderItems={dummyData.splicingLevelEst.value}
                    defaultVisible={true}
                >
                    <SplicingLevelSubtile data={dummyData.splicingLevelEst} onDimsChanged={this.props.onDimsChanged} />
                </CollapsibleSection>
            </CollapsibleTile>
        );
    }
};
