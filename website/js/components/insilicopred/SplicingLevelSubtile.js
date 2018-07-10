/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";
import classNames from 'classnames';
import {TabbedArea, TabPane} from "react-bootstrap";

import {capitalize, isNumeric} from "../../util";


// ================================================================
// === helper tables, functions
// ================================================================

// these impact tables are used by all the components (SplicingOverviewTable, SpliceSiteImpactTable, and
// DeNovoDonorPathogenicityTable) to color-code the pathogenicity probability

const splicingImpactFields = {
    fields: [
        { key: 'outside', label: 'Outside splice consensus region', prob: 0.02 },
        { key: 'improved', label: 'Improved', prob: 0.04 },
        { key: 'low', label: 'Weak/Null & Low', prob: 0.04 },
        { key: 'moderate', label: 'Moderate', prob: 0.34 },
        { key: 'high', label: 'High', prob: 0.97 }
    ],
    zScoreLabels: {
        donor: {
            outside: 'n/a',
            improved: 'Increased score',
            low: 'Z > 0',
            moderate: '-2.0 <= Z < 0.0',
            high: 'Z < -2.0'
        },
        acceptor: {
            outside: 'n/a',
            improved: 'Increased score',
            low: 'Z > 0.5',
            moderate: '-1.5 <= Z < 0.5',
            high: 'Z < -1.5'
        }
    }
};

const deNovoImpactFields = [
    {
        type: 'note',
        text: (<span>If this variant creates a de novo donor, the de novo donor creates an in-frame deletion that does <b><i>not</i></b> alter a key functional domain.</span>)
    },
    {
        type: 'value',
        impact: 'low',
        label: 'Innocuous IFD',
        zScoreLabel: 'n/a',
        prob: 0.02
    },
    {
        type: 'note',
        text: (<span>If this variant actually creates a de novo donor, the de novo donor either creates a frameshift or alters a key functional domain.</span>)
    },
    {
        type: 'value',
        impact: 'low',
        label: 'Weak/Null & Low',
        zScoreLabel: 'Z < -2',
        prob: 0.02
    },
    {
        type: 'value',
        impact: 'moderate',
        label: 'Moderate',
        zScoreLabel: '-2 >= Z > 0',
        prob: 0.30
    },
    {
        type: 'value',
        impact: 'high',
        label: 'Increased',
        zScoreLabel: 'Z > 0',
        prob: 0.64
    },
    {
        type: 'cond',
        options: [
            {
                /* special case #1: de novo donor */
                check: (data) => data.denovo.variant.zScore > data.denovo.closest.zScore,
                text: (<span>This variant was promoted one category because the Z score for the de novo splice site introduced by the variant is higher than the z score of the altered wild type donor sequence of this exon.</span>)
            },
            {
                check: () => true, // the default option, if reached
                text: (<span>If this variant creates a de novo donor, the de novo donor is not predicted to correct disruptions to the translation frame.</span>)
            }
        ]
    },
];

function getPathosLevel(prob, isDeNovo) {
    // FIXME: needs logic to disambiguate cases where the probability matches multiple rows

    if (!isDeNovo) {
        const elem = splicingImpactFields.fields.find(x => x.prob === prob);
        return elem ? elem.key : 'none';
    }
    else {
        const elem = deNovoImpactFields.find(x => x.prob === prob);
        return elem ? elem.impact : 'none';
    }
}

/**
 * Receives two strings, variant and ref, and compares each letter, bolding letters where they differ.
 * @param variant the sequence for the variant
 * @param ref the sequence for the reference
 * @returns {*} an array of React elements, either bolded if variant differs from ref or plain text if not.
 */
function boldedDiff(variant, ref) {
    // return variant;

    const matched = variant
        .split("")
        .slice(0, ref.length)
        .map((v, i) => ({match: v.toLowerCase() === ref[i].toLowerCase(), val: v }))
        .reduce((c, x) => {
            const last = c[c.length - 1];
            // collapse runs of matches or non-matches
            if (last && last.match === x.match) {
                last.val += x.val;
            }
            else {
                c.push(x);
            }

            return c;
        }, [])
        .map((x, i) => (x.match ? x.val : <span key={i} className="diff">{x.val}</span>));

    // return the array of elements with bolded differences, appending the remainder of variant if we didn't process it
    return (ref.length < variant.length) ? [matched, variant.slice(matched.length)] : matched;
}


// ================================================================
// === tab views
// ================================================================

class SplicingOverviewTable extends React.Component {
    render() {
        const {data, isDeNovo} = this.props;
        const variantRating = getPathosLevel(data.prior, isDeNovo);

        return (
            <table className="splicing-level-table" style={{marginTop: '1em'}}>
                <tbody>
                    <tr>
                        <td />
                        <td />
                        <td colSpan="2">MaxEntScan</td>
                        { isDeNovo && <td /> }
                    </tr>

                    <tr>
                        <td />
                        <td>Sequence</td>
                        <td>Raw Score</td>
                        <td>Z-score</td>
                        { isDeNovo && <td>Donor Position</td> }
                    </tr>

                    <tr>
                        <td>Wild Type</td>
                        <td className="sequence">{data.wild.sequence}</td>
                        <td>{data.wild.MES.toFixed(2)}</td>
                        <td>{data.wild.zScore.toFixed(2)}</td>
                        { isDeNovo &&  <td>{data.wild.donorPosition || 'n/a'}</td> }
                    </tr>

                    <tr>
                        <td>Variant</td>
                        <td className="sequence">{boldedDiff(data.variant.sequence, data.wild.sequence)}</td>
                        <td>{data.variant.MES.toFixed(2)}</td>
                        <td className={`pathos-prob-label-${variantRating}`} style={{fontWeight: 'bold'}}>
                        {data.variant.zScore.toFixed(2)}
                        </td>
                        { isDeNovo &&  <td>{data.variant.donorPosition || 'n/a'}</td> }
                    </tr>

                    {
                        isDeNovo && (
                            <tr>
                                <td>Reference Donor</td>
                                <td className="sequence">{data.closest.sequence}</td>
                                <td>{data.closest.MES.toFixed(2)}</td>
                                <td>{data.closest.zScore.toFixed(2)}</td>
                                <td>{data.closest.donorPosition}</td>
                            </tr>
                        )
                    }
                </tbody>
            </table>
        );
    }
}

class SpliceSiteImpactTable extends React.Component {
    render() {
        const {prior, type, variantZScore} = this.props;

        return (<div>
            <div style={{textAlign: 'center', fontWeight: 'bold', margin: '20px'}}>
                Probability of Pathogenicity<br />
                [Due to Splice {capitalize(type)} Damage]
            </div>

            <table className="splicing-level-table">
                <tr>
                    <td width="60%">Qualitative Category</td>
                    <td width="30%">Z-score Range</td>
                    <td width="20%">Probability</td>
                </tr>

                {
                    // FIXME: highlight row that matches the prior prob value
                    splicingImpactFields.fields.map(x =>
                        <tr key={x.key} className={prior === x.prob ? 'highlighted' : ''}>
                            <td className={`pathos-prob-label-${x.key}`}>{x.label}</td>
                            <td>{splicingImpactFields.zScoreLabels[type][x.key]}</td>
                            <td>{x.prob}</td>
                        </tr>
                    )
                }
            </table>

            {
                /* special case #2: native donor (when type == 'donor') */
                /* special case #3: native acceptor (when type == 'acceptor') */
                /* NOTE: they're combined because they use the same text, but might refactor later if they diverge */
                (prior === 0.34 && ((type === 'donor' && variantZScore < -2.0) || (type === 'acceptor' && variantZScore < -1.5)) && (
                    <div>
                    This variant was classified as Moderate because the Z score for the wildtype is already low and the change due to mutation is not considered large enough to warrant Severe.
                    </div>
                ))
            }
        </div>);
    }
}

class DeNovoDonorPathogenicityTable extends React.Component {
    render() {
        const {prior, data} = this.props;

        return (<div>
            <div style={{textAlign: 'center', fontWeight: 'bold', margin: '20px'}}>
                Probability of Pathogenicity<br />
                [Due to De Novo Donor]
            </div>

            <table className="splicing-level-table">
                <tr>
                    <td width="60%">Qualitative Category</td>
                    <td width="30%">Z-score Range</td>
                    <td width="20%">Probability</td>
                </tr>

                {
                    // FIXME: highlight row that matches the prior prob value
                    deNovoImpactFields.map((x, i) => {
                        switch(x.type) {
                            case 'value':
                                return (
                                    <tr key={i} className={prior === x.prob ? 'highlighted' : ''}>
                                        <td className={`pathos-prob-label-${x.impact}`}>{x.label}</td>
                                        <td>{x.zScoreLabel}</td>
                                        <td>{x.prob.toFixed(2)}</td>
                                    </tr>
                                );
                            case 'note':
                                return (
                                    <tr key={i}>
                                        <td colSpan={3} className="note-row">{x.text}</td>
                                    </tr>
                                );
                            case 'cond':
                                const found = x.options.find(x => x.check(data));
                                return (<tr key={i}>
                                    <td colSpan={3} className="note-row">
                                    { found ? found.text : 'n/a'}
                                    </td>
                                </tr>);
                            default:
                                return (<tr><td colSpan={3}>Unknown type: {x.type}</td></tr>);
                        }
                    })
                }
            </table>
        </div>);
    }
}


// ================================================================
// === splicing subtile implementation
// ================================================================

export default class SplicingLevelSubtile extends React.Component {
    constructor(props) {
        super(props);

        const maxProbPanel = this.getMaxProb();

        if (maxProbPanel) {
            this.state = {
                activePane: maxProbPanel.idx,
                initiallyActivePane: maxProbPanel.idx
            };
        }

        this.changePane = this.changePane.bind(this);
    }

    getMaxProb() {
        const data = this.props.data;
        const panes = [
            {idx: 0, prior: isNumeric(data.donor.prior) ? data.donor.prior : -9999,
                reason: 'splice donor damage' },
            {idx: 1, prior: isNumeric(data.denovo.prior) ? data.denovo.prior : -9999,
                reason: 'de novo donor splice site creation' },
            {idx: 2, prior: isNumeric(data.acceptor.prior) ? data.acceptor.prior : -9999,
                reason: 'splice acceptor damage' }
        ];
        const maxProb = Math.max(...panes.map(x => x.prior));

        if (isNumeric(maxProb) && maxProb !== -9999) {
            return panes.find(x => x.prior === maxProb);
        }
        else {
            return null;
        }
    }

    changePane(key) {
        this.setState({ activePane: key });
        this.props.onDimsChanged();
    }

    tabTitle(label, field, key) {
        return (
            <span className={classNames(
                "splicing-tab-header",
                key === this.state.initiallyActivePane && "highest-prob",
                !isNumeric(field) && "is-na",
            )}>
            {`${label} (${field})`}
            </span>
        );
    }

    render() {
        const {data} = this.props;

        const priorHeaders = {
            donor: isNumeric(data.donor.prior) ? data.donor.prior.toFixed(2) : 'N/A',
            denovo: isNumeric(data.denovo.prior) ? data.denovo.prior.toFixed(2) : 'N/A',
            acceptor: isNumeric(data.acceptor.prior) ? data.acceptor.prior.toFixed(2) : 'N/A'
        };

        const maxProbPanel = this.getMaxProb();

        if (!maxProbPanel) {
            // everything is N/A, so display a message corresponding to that
            return (
                <div className="subtile-container splicing-subtile" style={{padding: '0px'}}>
                    <div className="novalue-note">(No valid splicing-level estimation data was found.)</div>
                </div>
            );
        }

        return (
            <div className="subtile-container splicing-subtile" style={{padding: '0px'}}>
                <div className="preamble">
                The splicing-level estimation is due to <b>{maxProbPanel.reason}</b> which introduces a probability of pathogenicity of <b>{maxProbPanel.prior}</b>.
                </div>

                <TabbedArea activeKey={this.state.activePane} onSelect={this.changePane}>
                    <TabPane eventKey={0} tab={this.tabTitle("Donor Impact", priorHeaders.donor, 0)}>
                    {
                        isNumeric(data.donor.prior)
                            ? (
                                <div className="tab-body">
                                    <SplicingOverviewTable data={data.donor} />
                                    <SpliceSiteImpactTable prior={data.donor.prior} variantZScore={data.donor.variant.zScore} type='donor' />
                                </div>
                            )
                            : <div className="novalue-note">(replace with reasons why this might be N/A)</div>
                    }
                    </TabPane>

                    <TabPane eventKey={1} tab={this.tabTitle("De Novo Donor", priorHeaders.denovo, 1)}>
                    {
                        isNumeric(data.denovo.prior)
                            ? (
                                <div className="tab-body">
                                    <SplicingOverviewTable data={data.denovo} isDeNovo={true} />
                                    <DeNovoDonorPathogenicityTable prior={data.denovo.prior} data={data} />
                                </div>
                            )
                            : <div className="novalue-note">(replace with reasons why this might be N/A)</div>
                    }
                    </TabPane>

                    <TabPane eventKey={2} tab={this.tabTitle("Acceptor Impact", priorHeaders.acceptor, 2)}>
                    {
                        isNumeric(data.acceptor.prior)
                            ? (
                                <div className="tab-body">
                                    <SplicingOverviewTable data={data.acceptor} />
                                    <SpliceSiteImpactTable prior={data.acceptor.prior} variantZScore={data.acceptor.variant.zScore} type='acceptor' />
                                </div>
                            )
                            : <div className="novalue-note">(replace with reasons why this might be N/A)</div>
                    }
                    </TabPane>
                </TabbedArea>
            </div>
        );
    }
}
