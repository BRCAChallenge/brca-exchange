/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";
import {TabbedArea, TabPane} from "react-bootstrap";

function capitalize(w) {
    return w.charAt(0).toUpperCase() + w.substr(1);
}

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
        type: 'note',
        text: (<span>If this variant creates a de novo donor, the de novo donor is not predicted to correct disruptions to the translation frame.</span>)
    }
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
                        <td>{data.wild.MES}</td>
                        <td>{data.wild.zScore}</td>
                        { isDeNovo &&  <td>{data.wild.donorPosition || 'n/a'}</td> }
                    </tr>

                    <tr>
                        <td>Variant</td>
                        <td className="sequence">{data.variant.sequence}</td>
                        <td>{data.variant.MES}</td>
                        <td className={`pathos-prob-label-${variantRating}`} style={{fontWeight: 'bold'}}>{data.variant.zScore}</td>
                        { isDeNovo &&  <td>{data.variant.donorPosition || 'n/a'}</td> }
                    </tr>

                    {
                        isDeNovo && (
                            <tr>
                                <td>Reference Donor</td>
                                <td className="sequence">{data.closest.sequence}</td>
                                <td>{data.closest.MES}</td>
                                <td>{data.closest.zScore}</td>
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
        const {prior, type} = this.props;

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
        </div>);
    }
}

class DeNovoDonorPathogenicityTable extends React.Component {
    render() {
        const {prior} = this.props;

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
                    deNovoImpactFields.map((x, i) =>
                        x.type === 'value'
                            ? (
                                <tr key={i} className={prior === x.prob ? 'highlighted' : ''}>
                                    <td className={`pathos-prob-label-${x.impact}`}>{x.label}</td>
                                    <td>{x.zScoreLabel}</td>
                                    <td>{x.prob.toFixed(2)}</td>
                                </tr>
                            )
                            : (
                                <tr key={i}>
                                    <td colSpan={3} className="note-row">{x.text}</td>
                                </tr>
                            )
                    )
                }
            </table>
        </div>);
    }
}

function isNumeric(n) {
    return !isNaN(parseFloat(n)) && isFinite(n);
}

export default class SplicingLevelSubtile extends React.Component {
    constructor(props) {
        super(props);

        const maxProbPanel = this.getMaxProb();
        this.state = {
            activePane: maxProbPanel.idx
        };

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
        return panes.find(x => x.prior === maxProb);
    }

    changePane(key) {
        this.setState({ activePane: key });
        this.props.onDimsChanged();
    }

    render() {
        const {data} = this.props;

        const priorHeaders = {
            donor: isNumeric(data.donor.prior) ? data.donor.prior.toFixed(2) : 'N/A',
            denovo: isNumeric(data.denovo.prior) ? data.denovo.prior.toFixed(2) : 'N/A',
            acceptor: isNumeric(data.acceptor.prior) ? data.acceptor.prior.toFixed(2) : 'N/A'
        };

        const maxProbPanel = this.getMaxProb();

        return (
            <div className="subtile-container splicing-subtile" style={{padding: '0px'}}>
                <div className="preamble">
                The highest-probability explanation is <b>{maxProbPanel.reason}</b> with a probability of <b>{maxProbPanel.prior}</b>.
                The other tabs are included for more information.
                </div>

                <TabbedArea activeKey={this.state.activePane} onSelect={this.changePane}>
                    <TabPane eventKey={0} tab={`Donor Impact (${priorHeaders.donor})`}>
                    {
                        isNumeric(data.donor.prior)
                            ? (
                                <div className="tab-body">
                                    <SplicingOverviewTable data={data.donor} />
                                    <SpliceSiteImpactTable prior={data.donor.prior} type='donor' />
                                </div>
                            )
                            : <div className="novalue-note">(replace with reasons why this might be N/A)</div>
                    }
                    </TabPane>

                    <TabPane eventKey={1} tab={`De Novo Creation (${priorHeaders.denovo})`}>
                    {
                        isNumeric(data.denovo.prior)
                            ? (
                                <div className="tab-body">
                                    <SplicingOverviewTable data={data.denovo} isDeNovo={true} />
                                    <DeNovoDonorPathogenicityTable prior={data.denovo.prior} />
                                </div>
                            )
                            : <div className="novalue-note">(replace with reasons why this might be N/A)</div>
                    }
                    </TabPane>

                    <TabPane eventKey={2} tab={`Acceptor Impact (${priorHeaders.acceptor})`}>
                    {
                        isNumeric(data.acceptor.prior)
                            ? (
                                <div className="tab-body">
                                    <SplicingOverviewTable data={data.acceptor} />
                                    <SpliceSiteImpactTable prior={data.acceptor.prior} type='acceptor' />
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
