/*eslint-env browser */
/*global require: false, module */
'use strict';
import React from "react";
import util from "../../util";
import KeyInline from "../KeyInline";
import {Table} from "react-bootstrap";

class PathosProbScale extends React.Component {
    render() {
        const {value, brackets} = this.props;

        return (
            <svg xmlns="http://www.w3.org/2000/svg" width="100%" height="100" viewBox="0 0 500 60">
                <defs>
                    <linearGradient id="grad1" x1="0%" y1="0%" x2="100%" y2="0%">
                        <stop offset="0%"   style={{stopColor: '#0000FE', stopOpacity: 1}} />
                        <stop offset="100%" style={{stopColor: '#F60010', stopOpacity: 1}} />
                    </linearGradient>
                </defs>


                <g transform="translate(25 0)">
                    <rect x={0} y={0} width={450} height={40} strokeWidth={2} stroke="#ccc" fill="url(#grad1)" />

                    {
                        /* tickmarks */
                        Array.from({ length: brackets + 1 }).map((x, idx) => {
                            const val = (idx / brackets).toFixed(2);
                            return (
                                <text y={70} x={(450 * val)} alignmentBaseline="hanging"
                                    textAnchor={idx === 0 ? "start" : idx === brackets ? "end" : "middle"}>
                                    {val}
                                </text>
                            );
                        })
                    }

                    {/* indicators, if a value is given */}
                    { value !== undefined && (

                        <g>
                            <rect x={450 * value - 2} y={4} width={4} height={32} fill="black" />

                            <g transform={`translate(${450 * value} ${42})`}>
                            <polygon stroke="none" fill="black" points="0,0 5,12 -5,12" />
                            </g>
                        </g>
                    )}
                </g>
            </svg>
        );
    }
}


export default class InSilicoPredSubtile extends React.Component {
    render() {
        const {data} = this.props;

        const cols = data.context;
        // FIXME: the below will eventually come from an endpoint
        const summarizedData = data.context.reduce((c, a) => { c[a.prop] = a.value; return c; }, {});

        // for each panel, construct key-value pairs as a row of the table
        const submitterRows = cols.map(({prop, title, value, helpKey}) => {
            const isEmptyValue = util.isEmptyField(value);
            const rowItem = util.getFormattedFieldByProp(prop, summarizedData);

            return (
                <tr key={prop} className={ (isEmptyValue && this.props.hideEmptyItems) ? "variantfield-empty" : "" }>
                    {
                        helpKey
                            ? <KeyInline tableKey={title} onClick={(event) => this.props.showHelp(event, helpKey)}/>
                            : <td className='help-target'><span style={{fontWeight: 'bold'}}>{title}</span></td>
                    }
                    <td><span className="row-value">{rowItem}</span></td>
                </tr>
            );
        });

        return (
            <div>
                <div className="subtile-container" style={{padding: 0}}>
                    <PathosProbScale value={data.value} brackets={4} />
                </div>

                <div className="subtile-container">
                    Applicable Probability is from Splicing-level Estimation (De Novo Donor Splice Site Creation)
                </div>

                <Table style={{paddingBottom: 0}}>
                    <tbody>
                        {submitterRows}
                    </tbody>
                </Table>

                <div className="subtile-container">
                    <b>Credits:</b> Computational algorithm and display derived from the <a href="http://priors.hci.utah.edu/PRIORS/BRCA/">HCI Breast Cancer Genes Prior Probabilities</a> website.
                </div>
            </div>
        );
    }
}
