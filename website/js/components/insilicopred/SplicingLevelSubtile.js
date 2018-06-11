/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";
import {TabbedArea, TabPane} from "react-bootstrap";

function capitalize(w) {
    return w.charAt(0).toUpperCase() + w.substr(1);
}

class SplicingOverviewTable extends React.Component {
    render() {
        const {data, referenceDonor} = this.props;

        return (
            <table className="splicing-level-table" style={{marginTop: '1em'}}>
                <tbody>
                    <tr>
                        <td />
                        <td />
                        <td colSpan="2">MaxEntScan</td>
                        { referenceDonor && <td /> }
                    </tr>

                    <tr>
                        <td />
                        <td>Sequence</td>
                        <td>Raw Score</td>
                        <td>Z-score</td>
                        { referenceDonor && <td>Donor Position</td> }
                    </tr>

                    <tr>
                        <td>Wild Type</td>
                        <td className="sequence">{data.comparisons.wildType.sequence}</td>
                        <td>{data.comparisons.wildType.rawScore}</td>
                        <td>{data.comparisons.wildType.zScore}</td>
                        { referenceDonor &&  <td>{data.comparisons.wildType.donorPosition || 'n/a'}</td> }
                    </tr>

                    <tr>
                        <td>Variant</td>
                        <td className="sequence">{data.comparisons.variant.sequence}</td>
                        <td>{data.comparisons.variant.rawScore}</td>
                        <td>{data.comparisons.variant.zScore}</td>
                        { referenceDonor &&  <td>{data.comparisons.variant.donorPosition || 'n/a'}</td> }
                    </tr>

                    {
                        referenceDonor && (
                            <tr>
                                <td>Reference Donor</td>
                                <td className="sequence">{referenceDonor.sequence}</td>
                                <td>{referenceDonor.rawScore}</td>
                                <td>{referenceDonor.zScore}</td>
                                <td>{referenceDonor.donorPosition}</td>
                            </tr>
                        )
                    }
                </tbody>
            </table>
        );
    }
}

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
            moderate: '-2 <= Z < 0',
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

class SpliceSiteImpactTable extends React.Component {
    render() {
        const {type} = this.props;

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
                    splicingImpactFields.fields.map(x =>
                        <tr>
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

                <tr>
                    <td colSpan={3} className="note-row">
                    If this variant creates a de novo donor, the de novo donor creates an in-frame deletion that does <b><i>not</i></b> alter a key functional domain.
                    </td>
                </tr>

                <tr>
                    <td className="pathos-prob-label-low">Innocuous IFD</td>
                    <td>n/a</td>
                    <td>0.02</td>
                </tr>

                <tr>
                    <td colSpan={3} className="note-row">
                    If this variant actually creates a de novo donor, the de novo donor either creates a frameshift or alters a key functional domain.
                    </td>
                </tr>

                <tr>
                    <td className="pathos-prob-label-low">Weak/Null & Low</td>
                    <td>Z &lt; -2</td>
                    <td>0.02</td>
                </tr>

                <tr>
                    <td className="pathos-prob-label-moderate">Moderate</td>
                    <td>-2 &gt;= Z > 0</td>
                    <td>0.30</td>
                </tr>

                <tr>
                    <td className="pathos-prob-label-high">Increased</td>
                    <td>Z &gt; 0</td>
                    <td>0.64</td>
                </tr>

                <tr>
                    <td colSpan={3} className="note-row">
                    If this variant creates a de novo donor, the de novo donor is not predicted to correct disruptions to the translation frame.
                    </td>
                </tr>
            </table>
        </div>);
    }
}

export default class SplicingLevelSubtile extends React.Component {
    constructor(props) {
        super(props);

        this.state = {
            activePane: 1
        };

        this.changePane = this.changePane.bind(this);
    }

    changePane(key) {
        this.setState({ activePane: key });
        this.props.onDimsChanged();
    };

    render() {
        const {data} = this.props;

        const referenceDonor = {
            sequence: 'TGGgtaagg',
            rawScore: 9.08,
            zScore: 0.46,
            donorPosition: '547+1'
        };

        return (
            <div className="subtile-container splicing-subtile">
                <TabbedArea activeKey={this.state.activePane} onSelect={this.changePane}>
                    <TabPane eventKey={1} tab='Donor Impact'>
                        <SplicingOverviewTable data={data} />
                        <SpliceSiteImpactTable data={data} type='donor' />
                    </TabPane>

                    <TabPane eventKey={2} tab='De Novo Creation'>
                        <SplicingOverviewTable data={data} referenceDonor={referenceDonor} />
                        <DeNovoDonorPathogenicityTable data={data} type='donor' />
                    </TabPane>

                    <TabPane eventKey={3} tab='Acceptor Impact'>
                        <SplicingOverviewTable data={data} />
                        <SpliceSiteImpactTable data={data} type='acceptor' />
                    </TabPane>
                </TabbedArea>
            </div>
        );
    }
}
