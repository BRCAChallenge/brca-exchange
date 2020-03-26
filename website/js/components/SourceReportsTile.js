/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";
import {Panel} from 'react-bootstrap';
import util from '../util';
import slugify from '../slugify';
import {VariantSubmitter} from "./VariantSubmitter";

import GroupHelpButton from './GroupHelpButton';

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
        });
    }

    reportToggled(idx) {
        // return a new array in which only the selected element is toggled
        this.setState((pstate) => ({
            reportExpanded: pstate.reportExpanded.map((x, j) => (idx === j) ? !x : x)
        }));
    };

    render() {
        // put it in a temp b/c we're going to resort it
        let submissions = this.props.submissions;

        // get latest release id
        let latestReleaseID = 0;
        for (let i = 0; i < submissions.length; i++) {
            let releaseID = submissions[i].Data_Release.id;
            latestReleaseID = (releaseID > latestReleaseID) ? releaseID : latestReleaseID;
        }

        // filter out all old submissions
        let filteredSubmissions = submissions.filter(submission => submission.Data_Release.id === latestReleaseID);

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
                .map(({ title, prop, helpKey, noHelpLink }) => ({
                    title, prop, helpKey, noHelpLink, value: submissionData[prop]
                }));

            return (
                <VariantSubmitter
                    key={submissionData.id} idx={idx} submitter={submitterName} source={this.props.sourceName}
                    reportBinding={this.props.reportBinding} cols={formattedCols} data={submissionData}
                    hideEmptyItems={this.props.hideEmptyItems}
                    onReportToggled={this.reportToggled}
                    relayoutGrid={this.props.relayoutGrid}
                    showHelp={this.props.showHelp}
                    expanded={this.state.reportExpanded[idx]}
                    tooltips={this.props.tooltips}
                />
            );
        });

        // create the source panel itself now
        const groupTitle = `source-panel-${this.props.sourceName}`;

        return (
            <div key={`group_collection-${groupTitle}`} className={`variant-detail-group variant-submitter-group ${slugify(this.props.sourceName)}-submitter`}>
                <Panel
                    collapsible={true}
                    defaultExpanded={localStorage.getItem("collapse-group_" + groupTitle) !== "true"}
                >
                    <Panel.Heading>
                        <Panel.Title componentClass="h3">
                            <Panel.Toggle componentClass="a" className="title"
                                onClick={(event) => this.props.onChangeGroupVisibility(groupTitle, event)}
                            >
                                {this.props.groupTitle}
                            </Panel.Toggle>

                            <a title='collapse all reports'
                                className="toggle-subfields"
                                onClick={(event) => this.setAllReportExpansion(event, false)}
                                style={{cursor: 'pointer', marginRight: '10px'}}>
                                <i className="fa fa-angle-double-up" aria-hidden="true" />
                            </a>

                            <a title='expand all reports'
                                className="toggle-subfields"
                                onClick={(event) => this.setAllReportExpansion(event, true)}
                                style={{cursor: 'pointer'}}>
                                <i className="fa fa-angle-double-down" aria-hidden="true" />
                            </a>

                            {
                                this.props.helpSection &&
                                <GroupHelpButton group={this.props.helpSection}
                                    onClick={(event) => {
                                        this.props.showHelp(event, this.props.helpSection);
                                        return true;
                                    }}
                                />
                            }
                        </Panel.Title>
                    </Panel.Heading>
                    <Panel.Collapse
                        onEntered={this.props.relayoutGrid}
                        onExited={this.props.relayoutGrid}
                    >
                        <Panel.Body>
                        {submitters}
                        </Panel.Body>
                    </Panel.Collapse>
                </Panel>
            </div>
        );
    };
}
