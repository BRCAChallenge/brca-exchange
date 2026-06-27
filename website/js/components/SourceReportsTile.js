/*eslint-env browser */
'use strict';

import React from "react";
import { Card, Collapse } from 'react-bootstrap';
import util from '../util';
import slugify from '../slugify';
import { VariantSubmitter } from "./VariantSubmitter";
import GroupHelpButton from './GroupHelpButton';

export default class SourceReportsTile extends React.Component {
    constructor(props) {
        super(props);

        // seed group collapsed/expanded from localStorage (mirrors old defaultExpanded)
        const groupTitle = `source-panel-${props.sourceName}`;
        const defaultExpanded = localStorage.getItem("collapse-group_" + groupTitle) !== "true";

        this.state = {
            reportExpanded: this.defaultReportExpansions(),
            groupOpen: defaultExpanded
        };

        this.reportToggled = this.reportToggled.bind(this);
        this.setAllReportExpansion = this.setAllReportExpansion.bind(this);
        this.toggleGroupOpen = this.toggleGroupOpen.bind(this);
    }

    defaultReportExpansions() {
        // keep track of how many non-enigma/bic entries we've seen
        let nonEnigmaBics = 0;

        // temp because we may need to re-sort it
        let submissions = this.props.submissions;

        // sort the submissions if this source specifies a sort function
        if (this.props.reportBinding.sortBy) {
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
                nonEnigmaBics += 1;
            }

            // always collapse ENIGMA and BIC submissions.
            // show all items expanded if there are only a few of them.
            // otherwise, expand the first non-enigma/bic elem by default, but nothing else.
            return (!isEnigmaOrBic) && (this.props.submissions.length <= 3 || nonEnigmaBics === 1);
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

    toggleGroupOpen(event) {
        const groupTitle = `source-panel-${this.props.sourceName}`;
        if (this.props.onChangeGroupVisibility) {
            this.props.onChangeGroupVisibility(groupTitle, event);
        }
        this.setState(prev => ({ groupOpen: !prev.groupOpen }));
    }

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
                    key={submissionData.id}
                    idx={idx}
                    submitter={submitterName}
                    source={this.props.sourceName}
                    reportBinding={this.props.reportBinding}
                    cols={formattedCols}
                    data={submissionData}
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
            <div
                key={`group_collection-${groupTitle}`}
                className={`variant-detail-group variant-submitter-group ${slugify(this.props.sourceName)}-submitter`}
            >
                <Card className="mb-3 shadow">
                    <Card.Header as="div" className="d-flex align-items-center fw-bold">
                        <span
                            role="button"
                            className="title text-decoration-none flex-grow-1"
                            onClick={this.toggleGroupOpen}
                            aria-expanded={this.state.groupOpen}
                        >
                            {this.props.groupTitle}
                        </span>

                        <span
                            title='collapse all reports'
                            className="toggle-subfields"
                            onClick={(event) => this.setAllReportExpansion(event, false)}
                            style={{ cursor: 'pointer', marginRight: '10px' }}
                        >
                            <i className="fa fa-angle-double-up" aria-hidden="true" />
                        </span>

                        <span
                            title='expand all reports'
                            className="toggle-subfields"
                            onClick={(event) => this.setAllReportExpansion(event, true)}
                            style={{ cursor: 'pointer' }}
                        >
                            <i className="fa fa-angle-double-down" aria-hidden="true" />
                        </span>

                        {
                            this.props.helpSection &&
                            <GroupHelpButton
                                group={this.props.helpSection}
                                onClick={(event) => {
                                    this.props.showHelp(event, this.props.helpSection);
                                    return true;
                                }}
                            />
                        }
                    </Card.Header>

                    <Collapse
                        in={this.state.groupOpen}
                        onEntered={this.props.relayoutGrid}
                        onExited={this.props.relayoutGrid}
                    >
                        <div>
                            <Card.Body>
                                {submitters}
                            </Card.Body>
                        </div>
                    </Collapse>
                </Card>
            </div>
        );
    };
}

