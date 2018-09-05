/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";
import mapValues from 'lodash/mapValues';
import maxBy from 'lodash/maxBy';

import CollapsibleTile from "../collapsible/CollapsibleTile";
import CollapsibleSection from "../collapsible/CollapsibleSection";

import SplicingLevelSubtile from "./SplicingLevelSubtile";
import ProteinLevelSubtile from "./ProteinLevelSubtile";
import InSilicoPredSubtile from "./InSilicoPredSubtile";
import {isNumeric} from "../../util";

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
                MES: data.closestDonorAltMES,
                zScore: data.closestDonorAltZ,
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
        const noDataNode = (
            <CollapsibleTile allEmpty={true} {...this.props}>
                <div style={{padding: '10px'}}>
                No prior probability information exists for this variant.
                </div>
            </CollapsibleTile>
        );

        if (!this.props.priors) {
            return noDataNode;
        }

        // ensure any number-like fields are actually numbers (they're unfortunately stored as strings in the db to
        // handle '-', which is sometimes used to represent a missing value)
        const priorsData = mapValues(this.props.priors, (x) => isNumeric(x) ? parseFloat(x) : x);

        // ensure any NAs have been replaced with numeric values (i.e. -Infinity) before we do any comparisons
        const [proteinPrior, refDonorPrior, refAccPrior, deNovoDonorPrior] =
            [priorsData.proteinPrior, priorsData.refDonorPrior, priorsData.refAccPrior, priorsData.deNovoDonorPrior].map(numberify);

        // determine whether the probability of pathogenicity is the protein-level or splicing-level estimation
        // (de novo, acceptor, or donor).
        // whichever is highest is considered the deciding factor.
        const probAndReasons = [
            { prob: refDonorPrior, reason: 'Splicing-level Estimation (Donor Splice Site Impact)' },
            { prob: refAccPrior, reason: 'Splicing-level Estimation (Acceptor Splice Site Impact)' },
            { prob: deNovoDonorPrior, reason: 'Splicing-level Estimation (De Novo Donor Splice Site Creation)' },
            { prob: proteinPrior, reason: 'Protein-level Estimation' }
        ];
        const {prob: decidingProb, reason} = maxBy(probAndReasons, (x) => x.prob);

        // used by the splicing subtile to show the highest probability for that category
        const splicingPrior = Math.max(refDonorPrior, refAccPrior, deNovoDonorPrior);

        // restructure splicing-level props into a heirarchy so we don't have to pass each field separately
        const splicingData = extractSplicePayload(priorsData, true);

        // flag this entire tile as empty if none of the subtiles will contain valid data
        const allEmpty = [decidingProb, proteinPrior, splicingPrior].every(x => x === -Infinity);

        if (allEmpty) {
            // if we have no probabilities, the individual missing texts aren't accurate (e.g. protein probs
            // are just generally missing, not because the variant is intronic)
            return noDataNode;
        }

        return (
            <CollapsibleTile allEmpty={allEmpty} {...this.props}>
                <CollapsibleSection
                    fieldName={<span><i>In Silico</i> Prior Probability of Pathogenicity</span>}
                    extraHeaderItems={normalizedProb(decidingProb)}
                    defaultVisible={true}
                >
                    <InSilicoPredSubtile probability={decidingProb} reason={reason}
                        varLoc={priorsData.varLoc} />
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
