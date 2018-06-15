/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";
import classNames from "classnames";

const mapValToColor = (v) => {
    if (v < 0.25) {
        return "bracket_lt25";
    }
    else if (v < 0.5) {
        return "bracket_lt50";
    }
    else if (v < 0.75) {
        return "bracket_lt75";
    }
    else {
        return "bracket_lte100";
    }
};

const context = {
    missense: {
        keyDomainScores: [
            { name: 'C0', value: 0.03 },
            { name: 'C15-C25', value: 0.29 },
            { name: 'C35-C55', value: 0.66 },
            { name: 'C65', value: 0.81 },
            ],
        irrelevant: 0.02
    },
    silentVal: 0.02,
    nonsenseVal: 0.99
};

function labelClasses(v) {
    return classNames('label-field borderless ', mapValToColor(v));
}

function rowClasses(v, probability) {
    // FIXME: needs logic to disambiguate cases where the probability matches multiple rows
    return classNames('field-container', probability === v ? 'highlighted' : null);
}

export default class ProteinLevelSubtile extends React.Component {
    render() {
        const {probability} = this.props;

        return (
            <div className="subtile-container">
                <div className="protein-level-estimation" style={{border: 'solid 1px #eee', borderRight: 'none'}}>
                    <div className="field-container">
                        <div className="label-field borderless" />
                        <div className="value-field">
                            Probability
                        </div>
                    </div>

                    <div className="field-container">
                        <div className="label-field missense-header borderless">Missense (nonsynonymous)</div>
                        <div className="value-field" style={{border: 'none'}} />
                    </div>

                    <div className="nested">
                        <div className="field-container">
                            <div className="label-field" style={{border: 'none'}}>Inside key domains</div>
                        </div>

                        {
                            context.missense.keyDomainScores.map(x => (
                                <div className={rowClasses(x.value, probability)}>
                                    <div className={labelClasses(x.value)}>Align GV-GD score: {x.name}</div>
                                    <div className="value-field">{x.value}</div>
                                </div>
                            ))
                        }

                        <div className="field-container">
                            <div className="label-field" style={{border: 'none', marginTop: '0.5em'}}>Outside key domains</div>
                        </div>

                        <div className={rowClasses(context.missense.irrelevant, probability)}>
                            <div className={labelClasses(context.missense.irrelevant)}>
                                <i>Missense severity irrelevant</i>
                            </div>
                            <div className="value-field">{context.missense.irrelevant}</div>
                        </div>
                    </div>

                    <div className={rowClasses(context.silentVal, probability)} style={{marginTop: '1em'}}>
                        <div className={labelClasses(context.silentVal)}>
                            Silent (synonymous)
                        </div>
                        <div className="value-field">{context.silentVal}</div>
                    </div>

                    <div className={rowClasses(context.nonsenseVal, probability)} style={{marginTop: '1em'}}>
                        <div className={labelClasses(context.nonsenseVal)}>
                            Nonsense
                        </div>
                        <div className="value-field">{context.nonsenseVal}</div>
                    </div>
                </div>
            </div>
        );
    }
}
