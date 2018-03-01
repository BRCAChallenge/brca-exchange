/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";
import {Panel} from 'react-bootstrap';
import util from '../util';
import {VariantSubmitter} from "./VariantSubmitter";

export default class SourceReportsTile extends React.Component {
    constructor(props) {
        super(props);

        this.state = {
            reportExpanded: this.defaultReportExpansions()
        };

        this.reportToggled = this.reportToggled.bind(this);
        this.setAllReportExpansion = this.setAllReportExpansion.bind(this);
    }

    defaultReportExpansions() {
        // keep track of how many non-enigma/bic entries we've seen
        // (not a great solution b/c it introduces side effects into the map() below...)
        let nonEnigmaBics = 0;

        // put it in a temp b/c we may need to re-sort it
        // note that this is necessary because the order in which we process reports matters
        let submissions = this.props.submissions;

        // sort the submissions if this source specifies a sort function
        if (this.props.reportBinding.sortBy) {
            // (side note: we concat() to clone before sort()ing, because sort() mutates the array)
            submissions = submissions.concat().sort(this.props.reportBinding.sortBy);
        }

        return submissions.map((submissionData) => {
            const submitterName = util.getFormattedFieldByProp(this.props.reportBinding.submitter.prop, submissionData);

            const isEnigmaOrBic = (
                typeof submitterName === "string" &&
                (
                    submitterName.toLowerCase().indexOf("enigma") !== -1 ||
                    submitterName.toLowerCase().indexOf("(bic)") !== -1
                )
            );

            if (!isEnigmaOrBic) {
                // we only really care about the first, but this is the cleanest way to do this
                // with a single var
                nonEnigmaBics += 1;
            }

            // always collapse ENIGMA and BIC submissions.
            // show all items expanded if there are only a few of them.
            // otherwise, expand the first non-enigma/bic elem by default, but nothing else.
            return ( !isEnigmaOrBic ) && (this.props.submissions.length <= 3 || nonEnigmaBics === 1);
        });
    }

    setAllReportExpansion(e, newExpansion) {
        e.stopPropagation();

        this.setState({
            reportExpanded: Array.from({ length: this.props.submissions.length }, () => newExpansion)
        }, () => {
            // causes the parent to perform a (delayed) reflow
            this.props.onReportToggled();
        });
    }

    reportToggled(idx) {
        // return a new array in which only the selected element is toggled
        this.setState((pstate) => ({
            reportExpanded: pstate.reportExpanded.map((x, j) => (idx === j) ? !x : x)
        }), () => {
            // causes the parent to perform a (delayed) reflow
            this.props.onReportToggled();
        });
    };

    render() {
        // put it in a temp b/c we're going to resort it
        let submissions = this.props.submissions;

        // first, remove old versions of reports (they're only used for the diff)
        let seen = [];
        let filteredSubmissions = [];
        for (var i = 0; i < submissions.length; i++) {
            if (submissions[i].Source === "ClinVar") {
                var key = submissions[i].SCV_ClinVar;
            } else if (submissions[i].Source === "LOVD") {
                key = submissions[i].DBID_LOVD;
            }
            // reports are already sent to the ui sorted with the most recent first, so we can take the first one
            if (seen.indexOf(key) < 0) {
                filteredSubmissions.push(submissions[i]);
            }
            seen.push(key);
        }

        submissions = filteredSubmissions;

        // sort the submissions if this source specifies a sort function
        if (this.props.reportBinding.sortBy) {
            // (side note: we concat() to clone before sort()ing, because sort() mutates the array)
            submissions = submissions.concat().sort(this.props.reportBinding.sortBy);
        }

        // create a per-submitter collapsible subsection within this source panel
        const submitters = submissions.map((submissionData, idx) => {
            // extract header fields, e.g. the submitter name
            const submitterName = util.getFormattedFieldByProp(this.props.reportBinding.submitter.prop, submissionData);

            // extract fields we care about from the submission data
            const formattedCols = this.props.reportBinding.cols
                .map(({ title, prop, helpKey }) => ({
                    title, prop, helpKey, value: submissionData[prop]
                }));

            return (
                <VariantSubmitter
                    key={submissionData.id} idx={idx} submitter={submitterName} source={this.props.sourceName}
                    reportBinding={this.props.reportBinding} cols={formattedCols} data={submissionData}
                    hideEmptyItems={this.props.hideEmptyItems}
                    onReportToggled={this.reportToggled}
                    showHelp={this.props.showHelp}
                    expanded={this.state.reportExpanded[idx]}
                />
            );
        });

        // create the source panel itself now
        const groupTitle = `source-panel-${this.props.sourceName}`;
        const header = (
            <h3 style={{display: 'flex', flexDirection: 'row'}}>
                <a style={{flexGrow: 1}} href="#" onClick={(event) => this.props.onChangeGroupVisibility(groupTitle, event)}>
                    {this.props.groupTitle}
                </a>

                <a title='collapse all reports'
                    onClick={(event) => this.setAllReportExpansion(event, false)}
                    style={{cursor: 'pointer', marginRight: '10px'}}>
                    <i className="fa fa-angle-double-up" aria-hidden="true" />
                </a>

                <a title='expand all reports'
                    onClick={(event) => this.setAllReportExpansion(event, true)}
                    style={{cursor: 'pointer'}}>
                    <i className="fa fa-angle-double-down" aria-hidden="true" />
                </a>
            </h3>
        );

        return (
            <div key={`group_collection-${groupTitle}`} className="variant-detail-group variant-submitter-group">
                <Panel
                    header={header}
                    collapsable={true}
                    defaultExpanded={localStorage.getItem("collapse-group_" + groupTitle) !== "true"}>
                    {submitters}
                </Panel>
            </div>
        );
    };
}
