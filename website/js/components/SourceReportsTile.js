/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";
import {Panel} from 'react-bootstrap';
import util from '../util';
import {reportSourceFieldMapping} from "../VariantTable";
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
        // get metadata about source to determine if we need to sort, which prop is the submitter name, etc.
        const sourceMeta = reportSourceFieldMapping[this.props.sourceName];

        // keep track of how many non-enigma/bic entries we've seen
        // (not a great solution b/c it introduces side effects into the map() below...)
        let nonEnigmaBics = 0;

        // put it in a temp b/c we may need to re-sort it
        // note that this is necessary because the order in which we process reports matters
        let submissions = this.props.submissions;

        // sort the submissions if this source specifies a sort function
        if (sourceMeta.sortBy) {
            // (side note: we concat() to clone before sort()ing, because sort() mutates the array)
            submissions = submissions.concat().sort(sourceMeta.sortBy);
        }

        return submissions.map((submissionData) => {
            const submitterName = util.getFormattedFieldByProp(sourceMeta.submitter.prop, submissionData);

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
        // look up data for formatting this specific source, e.g. what fields to include, its name, etc.
        const sourceMeta = reportSourceFieldMapping[this.props.sourceName];

        // put it in a temp b/c we're going to resort it
        let submissions = this.props.submissions;

        // sort the submissions if this source specifies a sort function
        if (sourceMeta.sortBy) {
            // (side note: we concat() to clone before sort()ing, because sort() mutates the array)
            submissions = submissions.concat().sort(sourceMeta.sortBy);
        }

        // create a per-submitter collapsible subsection within this source panel
        const submitters = submissions.map((submissionData, idx) => {
            // extract header fields, e.g. the submitter name
            const submitterName = util.getFormattedFieldByProp(sourceMeta.submitter.prop, submissionData);

            // extract fields we care about from the submission data
            const formattedCols = sourceMeta.cols
                .map(({ title, prop }) => ({
                    title, prop, value: submissionData[prop]
                        ? submissionData[prop].toString()
                        // FIXME: maybe we should just ignore requested fields that aren't in the data payload
                        : `${prop}???`
                }));

            return (
                <VariantSubmitter
                    key={submissionData.id} idx={idx} submitter={submitterName} source={this.props.sourceName}
                    meta={sourceMeta} cols={formattedCols} data={submissionData}
                    hideEmptyItems={this.state.hideEmptyItems}
                    onReportToggled={this.reportToggled}
                    showHelp={this.showHelp}
                    expanded={this.state.reportExpanded[idx]}
                />
            );
        });

        // create the source panel itself now
        const groupTitle = `source-panel-${this.props.sourceName}`;
        const header = (
            <h3 style={{display: 'flex', flexDirection: 'row'}}>
                <a style={{flexGrow: 1}} href="#" onClick={(event) => this.props.onChangeGroupVisibility(groupTitle, event)}>
                    {`Clinical Significance (${sourceMeta.displayName})`}
                </a>

                <a title='collapse all reports'
                    onClick={(event) => this.setAllReportExpansion(event, false)}
                    style={{cursor: 'pointer', color: (this.state.reportExpanded.every(x => !x) ? 'black' : 'gray'), marginRight: '10px'}}>
                    <i className="fa fa-angle-double-up" aria-hidden="true" />
                </a>

                <a title='expand all reports'
                    onClick={(event) => this.setAllReportExpansion(event, true)}
                    style={{cursor: 'pointer', color: (this.state.reportExpanded.every(x => x) ? 'black' : 'gray')}}>
                    <i className="fa fa-angle-double-down" aria-hidden="true" />
                </a>
            </h3>
        );
        const allEmpty = false; // FIXME: actually check if we're all empty or no

        return (
            <div key={`group_collection-${groupTitle}`} className="variant-detail-group variant-submitter-group">
                <div className={ allEmpty && this.state.hideEmptyItems ? "group-empty" : "" }>
                    <Panel
                        header={header}
                        collapsable={true}
                        defaultExpanded={localStorage.getItem("collapse-group_" + groupTitle) !== "true"}>
                        {submitters}
                    </Panel>
                </div>
            </div>
        );
    };
}
