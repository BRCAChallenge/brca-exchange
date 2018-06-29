/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";
import {CollapsableMixin, Table} from "react-bootstrap";
import classNames from 'classnames';
import util from '../util';
import KeyInline from './KeyInline';

import slugify from '../slugify';

// from https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/#revstat_def
const marksToReviewStatuses = {
    4: ['practice_guideline'],
    3: ['reviewed_by_expert_panel'],
    2: [],
    1: ['criteria_provided,_single_submitter'],
    0: ['no_assertion_criteria_provided', 'no_assertion_provided'],
};

function computeReviewStatusScore(status) {
    // determine how many marks correspond to this review status
    for (let marksStatus of Object.entries(marksToReviewStatuses)) {
        if (marksStatus[1].includes(status)) {
            return marksStatus[0];
        }
    }

    return 0;
}

function getMarksForReviewStatus(status) {
    const marks = computeReviewStatusScore(status);
    const icon = 'star';

    // returns an array of length 4 (the maximum number of marks) where
    // the number of marks corresponding to this status are highlighted
    return (
        <span title={`${marks} ${icon} out of 4`} style={{display: 'inline-block', marginRight: '10px'}}>
        {
            Array.from({ length: 4 })
                .map((_, x) => (
                    <i key={x} className={`fa fa-${icon}`}
                        style={{ color: (x < marks) ? '#ffe86d' : '#cbcbcb', marginRight: '2px' }}
                        aria-hidden="true"
                    />
                ))
        }
        </span>
    );
}

const VariantSubmitter = React.createClass({
    mixins: [CollapsableMixin],

    getCollapsableDOMNode: function() {
        return this.refs.panel.getDOMNode();
    },

    getCollapsableDimensionValue: function() {
        return this.refs.panel.getDOMNode().scrollHeight;
    },

    onHandleToggle: function(e) {
        e.preventDefault();

        // ask our parent to toggle us
        this.props.onReportToggled(this.props.idx);
    },

    // source-specific headers
    generateHeader: function(source, submitter, data) {
        let extraHeaderItems = null;

        if (source === 'ClinVar') {
            const significance = util.getFormattedFieldByProp("Clinical_Significance_ClinVar", data)
                .replace(/(variant of unknown significance|uncertain significance)/i, 'VUS');
            const dateUpdated = util.getFormattedFieldByProp("Date_Last_Updated_ClinVar", data);

            extraHeaderItems = (
                <div className="clinvar-extras">
                    <span className="clin-sig">{util.sentenceCase(significance)}</span>
                    <span className="date-updated">{` (${dateUpdated})`}</span>
                </div>
            );
        }
        else if (source === 'LOVD') {
            const variantEffect = util.getFormattedFieldByProp("Variant_effect_LOVD", data);

            extraHeaderItems = (
                <div style={{whiteSpace: 'nowrap', overflow: 'hidden'}}>
                {`Variant Effect: ${variantEffect}`}
                </div>
            );
        }

        return (
            <div className={`submitter-header ${this.props.expanded ? 'expanded' : ''}`} onClick={this.onHandleToggle}>
                <div className="submitter-cell submitter-label">
                    {
                        // for ClinVar reports, display stars to the left of the submitter label
                        (this.props.source === 'ClinVar')
                            ? getMarksForReviewStatus(data['Review_Status_ClinVar'])
                            : null
                    }

                    {
                        this.props.expanded
                            ? <i className="fa fa-caret-down" aria-hidden="true" />
                            : <i className="fa fa-caret-right" aria-hidden="true" />
                    }
                    &nbsp;
                    <span>{this.props.reportBinding.submitter.title}</span>
                </div>

                <div
                    className={`submitter-cell submitter-name ${!this.props.expanded ? 'collapsed' : ''}`}
                    style={{textAlign: 'left', flex: '1 1 auto'}}
                >
                { util.abbreviatedSubmitter(submitter) }
                </div>

                <div className="submitter-cell optional" style={{textAlign: 'left', flex: '0 1 auto'}}>
                {
                    // remaining header elements depend on the source
                    extraHeaderItems
                }
                </div>
            </div>
        );
    },

    render: function() {
        let styles = this.getCollapsableClassSet();
        const {submitter, cols, data} = this.props;

        // for each panel, construct key-value pairs as a row of the table
        const submitterRows = cols.map(({prop, title, value, noHelpLink}) => {
            const isEmptyValue = util.isEmptyField(value);
            const rowItem = util.getFormattedFieldByProp(prop, data);

            return (
                <tr key={prop} className={ (isEmptyValue && this.props.hideEmptyItems) ? "variantfield-empty" : "" }>
                    {
                        !noHelpLink
                            ? (
                                <KeyInline tableKey={title}
                                    tooltip={this.props.tooltips && this.props.tooltips["report-" + slugify(prop)]}
                                    onClick={(event) => this.props.showHelp(event, "report-" + prop)}/>
                            )
                            : <td className='help-target'><span style={{fontWeight: 'bold'}}>{title}</span></td>
                    }
                    <td><span className="row-value">{rowItem}</span></td>
                </tr>
            );
        });

        return (
            <div>
                <div style={{marginBottom: 0, borderTop: 'solid 2px #ccc'}}>
                {
                    // remaining header elements depend on the source
                    this.generateHeader(this.props.source, submitter, data)
                }
                </div>

                <div ref='panel' className={classNames(styles)}>
                    <Table key={`submitter-name-${submitter}`}>
                        <tbody>
                        { submitterRows }
                        </tbody>
                    </Table>
                </div>
            </div>
        );
    }
});

module.exports = {
    VariantSubmitter,
    computeReviewStatusScore
};
