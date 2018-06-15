/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";

import CollapsibleTile from "../collapsible/CollapsibleTile";
import CollapsibleSection from "../collapsible/CollapsibleSection";

import SplicingLevelSubtile from "./SplicingLevelSubtile";
import ProteinLevelSubtile from "./ProteinLevelSubtile";
import InSilicoPredSubtile from "./InSilicoPredSubtile";

const mockData = {
    varLoc: 'Exon',
    applicablePrior: 0.03, // derived from max(proteinPrior, max(refDonorPrior, deNovoDonorPrior, refAccPrior))
    proteinPrior: 0.02,

    refDonorPrior: 0.04,
    deNovoDonorPrior: 0.3,
    refAccPrior: 'n/a',

    // due to donor damage
    refRefDonorSeq: 'TTGgtaaaa',
    refRefDonorMES: 3.23,
    refRefDonorZ: -3.293,
    altRefDonorSeq: 'TTGgaaaaa',
    altRefDonorMES: -4.95,
    altRefDonorZ: -5.6,

    // due to de novo donor
    refDeNovoDonorSeq: 'ACTGTGAGA',
    refDeNovoDonorMES: -0.91,
    refDeNovoDonorZ: -3.80,
    altDeNovoDonorSeq: 'ACGgtgaga',
    altDeNovoDonorMES: 5.64,
    altDeNovoDonorZ: 0.4,
    deNovoDonorTranscriptSplicePos: '64',
    deNovoDonorGenomicSplicePos: '64',

    closestDonorRefSeq: 'TGGgtaagg',
    closestDonorRefMES: 9.08,
    closestDonorRefZ: 0.46,
    closestDonorAltSeq: 'TAGgtattg',
    closestDonorAltMES: 5.46,
    closestDonorAltZ: -0.99,
    closestDonorTranscriptSplicePos: '67 + 1',
    closestDonorGenomicSplicePos: '67 + 1',

    // due to acceptor damage
    refRefAccSeq: 'ttttcccttgtattttacagATG',
    refRefAccMES: 11.68,
    refRefAccZ: 1.51,
    altRefAccSeq: 'ttttcccttgtattgtacagATG',
    altRefAccMES: 8.91,
    altRefAccZ: 0.38,
};

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

export default class SilicoPredTile extends React.Component {
    constructor(props) {
        super(props);
    }

    render() {
        // determine whether the probability of pathogenicity is the protein-level or splicing-level
        // estimation (de novo, acceptor, or donor).
        // whichever is highest is considered the deciding factor.
        const splicingPrior = Math.max(mockData.refDonorPrior, mockData.deNovoDonorPrior);
        const decidingProb = Math.max(mockData.proteinPrior, splicingPrior);
        const reason = mockData.proteinPrior > splicingPrior
            ? 'Protein-level Estimation'
            : 'Splicing-level Estimation ' + (
                mockData.deNovoDonorPrior > mockData.refDonorPrior
                    ? '(De Novo Donor Splice Site Creation)'
                    : `(${mockData.refDonorPrior > mockData.refAccPrior ? 'Donor' : 'Acceptor'} Splice Site Impact)`
            );

        /*
        notes from slide: "Variant Type: varType is NOT correct; only indicates insertion, deletion, delins, substitution.
        Can be assigned with proteinPrior as long as location relative to CI domain is accounted for (i.e. check that 0.02
        is not assigned bc of outside CI domain to assign a variant as silent)."
         */
        const varType = '???'; // FIXME: how do we determine the variant type?

        // restructure splicing-level props into a heirarchy so we don't have to pass them separately
        const splicingData = extractSplicePayload(mockData, true);

        return (
            <CollapsibleTile {...this.props}>
                <CollapsibleSection
                    fieldName="In Silico Probability of Pathogenicity"
                    extraHeaderItems={decidingProb}
                    defaultVisible={true}
                >
                    <InSilicoPredSubtile probability={decidingProb} reason={reason}
                        varLoc={mockData.varLoc} varType={varType} />
                </CollapsibleSection>

                <CollapsibleSection
                    fieldName="Protein-level Estimation"
                    extraHeaderItems={mockData.proteinPrior}
                    defaultVisible={false}
                >
                    <ProteinLevelSubtile probability={mockData.proteinPrior} />
                </CollapsibleSection>

                <CollapsibleSection
                    fieldName="Splicing-level Estimation"
                    extraHeaderItems={splicingPrior}
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
