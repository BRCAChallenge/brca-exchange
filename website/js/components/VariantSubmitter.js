/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";
import {CollapsableMixin, Table} from "react-bootstrap";
import classNames from 'classnames';
import util from '../util';
import KeyInline from './KeyInline';

function sentenceCase(str) {
    return str.replace(/\b\S/g, (t) => t.toUpperCase() );
}

// from https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/#revstat_def
const starsToRevstatuses = {
    4: ['practice_guideline'],
    3: ['reviewed_by_expert_panel'],
    2: [],
    1: ['criteria_provided,_single_submitter'],
    0: ['no_assertion_criteria_provided', 'no_assertion_provided'],
};

function getStarsForReviewStatus(status) {
    let stars = 0;

    // determine how many stars correspond to this review status
    for (let starStatus of Object.entries(starsToRevstatuses)) {
        if (starStatus[1].includes(status)) {
            stars = starStatus[0];
            break;
        }
    }

    // returns an array of length 4 (the maximum number of stars) where
    // the number of stars corresponding to this status are highlighted
    return (<span title={`${stars} stars out of 4`} style={{display: 'inline-block', marginRight: '10px'}}>
    {
        Array.from({ length: 4 })
            .map((_, x) => (
                <i key={x} className="fa fa-star"
                    style={{ color: (x < stars) ? '#ffe86d' : '#cbcbcb' }}
                    aria-hidden="true"
                />
            ))
    }
    </span>);
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
        this.setState({expanded: !this.state.expanded});
    },

    // FIXME: copied from index.js; find a way to use this without aliasing
    truncateData: function(field) {
        const fieldsToTruncate = ["Genomic_Coordinate_hg38", "Genomic_Coordinate_hg37", "Genomic_Coordinate_hg36"];
        return fieldsToTruncate.indexOf(field) > -1;
    },

    // source-specific headers
    generateHeader: function(source, submitter, data) {
        let extraHeaderItems = null;

        if (source === 'ClinVar') {
            const significance = util.getFormattedFieldByProp("Clinical_Significance_ClinVar", data);
            const dateUpdated = util.getFormattedFieldByProp("Date_Last_Updated_ClinVar", data);

            extraHeaderItems = (
                <div>
                {`${sentenceCase(significance)} (${dateUpdated})`}
                </div>
            );
        }
        else if (source === 'LOVD') {
            const variantEffect = util.getFormattedFieldByProp("Variant_effect_LOVD", data);

            extraHeaderItems = (
                <div>
                {`Variant Effect: ${variantEffect}`}
                </div>
            );
        }

        return (
            <div className={`submitter-header ${this.state.expanded ? 'expanded' : ''}`} onClick={this.onHandleToggle}>
                <div className="submitter-cell submitter-label">
                    {
                        // for ClinVar reports, display stars to the left of the submitter label
                        (this.props.source === 'ClinVar')
                            ? getStarsForReviewStatus(data['Review_Status_ClinVar'])
                            : null
                    }

                    {
                        this.state.expanded
                            ? <i className="fa fa-caret-down" aria-hidden="true" />
                            : <i className="fa fa-caret-right" aria-hidden="true" />
                    }
                    &nbsp;
                    <span>{this.props.meta.submitter.title}</span>
                </div>

                <div
                    className={`submitter-cell submitter-name ${!this.state.expanded ? 'collapsed' : ''}`}
                    style={{textAlign: 'left', flex: '1 1 auto'}}
                >
                { submitter }
                </div>

                <div className="submitter-cell" style={{textAlign: 'left', flex: '0 1 auto', whiteSpace: 'nowrap'}}>
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

        // ----

        const {submitter, cols, data} = this.props;

        // for each panel, construct key-value pairs as a row of the table
        const submitterRows = cols.map(({prop, title, value}) => {
            const isEmptyValue = util.isEmptyField(value);
            const rowItem = util.getFormattedFieldByProp(prop, data);

            return (
                <tr key={prop} className={ (isEmptyValue && this.props.hideEmptyItems) ? "variantfield-empty" : "" }>
                    <KeyInline tableKey={title} onClick={(event) => this.props.showHelp(event, title)} />
                    <td colSpan={2} ><span className={ this.truncateData(prop) ? "row-value-truncated" : "row-value" }>{rowItem}</span></td>
                    <td>&nbsp;</td>
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

export default VariantSubmitter;
