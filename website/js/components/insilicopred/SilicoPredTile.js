/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";

import CollapsibleTile from "../collapsible/CollapsibleTile";
import CollapsibleSection from "../collapsible/CollapsibleSection";

import SplicingLevelSubtile from "./SplicingLevelSubtile";
import ProteinLevelSubtile from "./ProteinLevelSubtile";
import InSilicoPredSubtile from "./InSilicoPredSubtile";
import {isNumeric} from "../../util";

const examples = JSON.parse(require('raw-loader!./insilico_mock.min.json'));
import {Panel} from "react-bootstrap";

function extractSplicePayload(data, useTranscriptSplicePos) {
    return {
        donor: {
            prior: data.refDonorPrior,
            wild: {
                sequence: data.refRefDonorSeq,
                MES: data.refRefDonorMES,
                zScore: data.refRefDonorZ
            },
            variant: {
                sequence: data.altRefDonorSeq,
                MES: data.altRefDonorMES,
                zScore: data.altRefDonorZ
            }
        },
        denovo: {
            prior: data.deNovoDonorPrior,
            wild: {
                sequence: data.refDeNovoDonorSeq,
                MES: data.refDeNovoDonorMES,
                zScore: data.refDeNovoDonorZ,
                donorPosition: useTranscriptSplicePos ? data.deNovoDonorTranscriptSplicePos : data.deNovoDonorGenomicSplicePos
            },
            variant: {
                sequence: data.altDeNovoDonorSeq,
                MES: data.altDeNovoDonorMES,
                zScore: data.altDeNovoDonorZ,
                donorPosition: useTranscriptSplicePos ? data.deNovoDonorTranscriptSplicePos : data.deNovoDonorGenomicSplicePos
            },
            closest: {
                sequence: data.closestDonorRefSeq,
                MES: data.closestDonorRefMES,
                zScore: data.closestDonorRefZ,
                donorPosition: useTranscriptSplicePos ? data.closestDonorTranscriptSplicePos : data.closestDonorGenomicSplicePos
            },
            closestAlt: {
                sequence: data.closestDonorAltSeq,
                MES: data.closestDonorAltMES, // where is this?
                zScore: data.closestDonorAltZ, // where is this?
                donorPosition: useTranscriptSplicePos ? data.closestDonorTranscriptSplicePos : data.closestDonorGenomicSplicePos
            }
        },
        acceptor: {
            prior: data.refAccPrior,
            wild: {
                sequence: data.refRefAccSeq,
                MES: data.refRefAccMES,
                zScore: data.refRefAccZ
            },
            variant: {
                sequence: data.altRefAccSeq,
                MES: data.altRefAccMES,
                zScore: data.altRefAccZ
            }
        }
    };
}

function numberify(x) {
    return isNumeric(x) ? x : -Infinity;
}

function normalizedProb(x) {
    return x !== -Infinity ? x : 'N/A';
}

export default class SilicoPredTile extends React.Component {
    constructor(props) {
        super(props);
    }

    render() {
        // extract the corresponding entry from the mock data, if it exists
        const mockData = examples[this.props.Genomic_Coordinate_hg38];

        if (!mockData) {
            const groupVisID = `group-panel-${this.props.groupTitle}`;
            const allEmpty = false;
            return (
                <div key={`group_collection-${groupVisID}`} className={ allEmpty && this.props.hideEmptyItems ? "group-empty variant-detail-group" : "variant-detail-group" }>
                    <Panel header={<h3><i>In Silico</i> Prediction</h3>}>
                        <div style={{padding: '10px'}}>
                        no mock data found for variant
                        </div>
                    </Panel>
                </div>
            );
        }

        // determine whether the probability of pathogenicity is the protein-level or splicing-level
        // estimation (de novo, acceptor, or donor).
        // whichever is highest is considered the deciding factor.

        // ensure any NAs have been replaced with numeric values (i.e. -Infinity) before we do any comparisons
        const [proteinPrior, refDonorPrior, refAccPrior, deNovoDonorPrior] =
            [mockData.proteinPrior, mockData.refDonorPrior, mockData.refAccPrior, mockData.deNovoDonorPrior].map(numberify);

        const splicingPrior = Math.max(refDonorPrior, refAccPrior, deNovoDonorPrior);
        const decidingProb = Math.max(proteinPrior, splicingPrior);
        const reason = proteinPrior > splicingPrior
            ? 'Protein-level Estimation'
            : 'Splicing-level Estimation ' + (
                deNovoDonorPrior > refDonorPrior
                    ? '(De Novo Donor Splice Site Creation)'
                    : `(${refDonorPrior > refAccPrior ? 'Donor' : 'Acceptor'} Splice Site Impact)`
            );

        /*
        notes from slide: "Variant Type: varType is NOT correct; only indicates insertion, deletion, delins, substitution.
        Can be assigned with proteinPrior as long as location relative to CI domain is accounted for (i.e. check that 0.02
        is not assigned bc of outside CI domain to assign a variant as silent)."
         */
        const varType = '???'; // FIXME: how do we determine the variant type?

        // restructure splicing-level props into a heirarchy so we don't have to pass them separately
        const splicingData = extractSplicePayload(mockData, true);

        // flag this entire tile as empty if none of the subtiles will contain valid data
        const allEmpty = [decidingProb, proteinPrior, splicingPrior].every(x => x === -Infinity);

        return (
            <CollapsibleTile allEmpty={allEmpty} {...this.props}>
                <CollapsibleSection
                    fieldName={<span><i>In Silico</i> Probability of Pathogenicity</span>}
                    extraHeaderItems={normalizedProb(decidingProb)}
                    defaultVisible={true}
                >
                    <InSilicoPredSubtile probability={decidingProb} reason={reason}
                        varLoc={mockData.varLoc} varType={varType} />
                </CollapsibleSection>

                <CollapsibleSection
                    fieldName="Protein-level Estimation"
                    extraHeaderItems={normalizedProb(proteinPrior)}
                    defaultVisible={false}
                >
                    <ProteinLevelSubtile probability={proteinPrior} />
                </CollapsibleSection>

                <CollapsibleSection
                    fieldName="Splicing-level Estimation"
                    extraHeaderItems={normalizedProb(splicingPrior)}
                    defaultVisible={false}
                >
                    <SplicingLevelSubtile
                        data={splicingData}
                        onDimsChanged={this.props.onDimsChanged}
                    />
                </CollapsibleSection>
            </CollapsibleTile>
        );
    }
};
