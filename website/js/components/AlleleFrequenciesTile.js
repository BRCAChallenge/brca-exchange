/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";
import {Panel, Table} from 'react-bootstrap';
import util from '../util';
const _ = require('underscore');
const KeyInline = require('./KeyInline');


export default class AlleleFrequenciesTile extends React.Component {
    constructor(props) {
        super(props);

        this.state = {
            // reportExpanded: this.defaultReportExpansions()
        };

        // this.reportToggled = this.reportToggled.bind(this);
        // this.setAllReportExpansion = this.setAllReportExpansion.bind(this);
    }

    getRows(source, data, variant) {
        return _.map(data, (rowDescriptor) => {
            let {prop, title} = rowDescriptor;
            let rowItem;

            if (variant[prop] !== null) {
                rowItem = util.getFormattedFieldByProp(prop, variant);
            }

            let isEmptyValue = util.isEmptyField(variant[prop]);

            if (isEmptyValue) {
                // rowsEmpty += 1;
                rowItem = '-';
            }

            // totalRowsEmpty += rowsEmpty;
            return (
                <tr key={prop} className={ (isEmptyValue && this.state.hideEmptyItems) ? "variantfield-empty" : "" }>
                    { rowDescriptor.tableKey !== false &&
                        <KeyInline tableKey={title} onClick={(event) => this.props.showHelp(event, title)}/>
                    }
                    <td colSpan={rowDescriptor.tableKey === false ? 2 : null} ><span className={"row-value" }>{rowItem}</span></td>
                </tr>
            );
        });
    }
    // defaultReportExpansions() {
        // // keep track of how many non-enigma/bic entries we've seen
        // // (not a great solution b/c it introduces side effects into the map() below...)
        // let nonEnigmaBics = 0;

        // // put it in a temp b/c we may need to re-sort it
        // // note that this is necessary because the order in which we process reports matters
        // let submissions = this.props.submissions;

        // // sort the submissions if this source specifies a sort function
        // if (this.props.reportBinding.sortBy) {
        //     // (side note: we concat() to clone before sort()ing, because sort() mutates the array)
        //     submissions = submissions.concat().sort(this.props.reportBinding.sortBy);
        // }

        // return submissions.map((submissionData) => {
        //     const submitterName = util.getFormattedFieldByProp(this.props.reportBinding.submitter.prop, submissionData);

        //     const isEnigmaOrBic = (
        //         typeof submitterName === "string" &&
        //         (
        //             submitterName.toLowerCase().indexOf("enigma") !== -1 ||
        //             submitterName.toLowerCase().indexOf("(bic)") !== -1
        //         )
        //     );

        //     if (!isEnigmaOrBic) {
        //         // we only really care about the first, but this is the cleanest way to do this
        //         // with a single var
        //         nonEnigmaBics += 1;
        //     }

        //     // always collapse ENIGMA and BIC submissions.
        //     // show all items expanded if there are only a few of them.
        //     // otherwise, expand the first non-enigma/bic elem by default, but nothing else.
        //     return ( !isEnigmaOrBic ) && (this.props.submissions.length <= 3 || nonEnigmaBics === 1);
        // });
    // }

    // setAllReportExpansion(e, newExpansion) {
    //     e.stopPropagation();

    //     this.setState({
    //         reportExpanded: Array.from({ length: this.props.submissions.length }, () => newExpansion)
    //     }, () => {
    //         // causes the parent to perform a (delayed) reflow
    //         this.props.onReportToggled();
    //     });
    // }

    // reportToggled(idx) {
    //     // return a new array in which only the selected element is toggled
    //     this.setState((pstate) => ({
    //         reportExpanded: pstate.reportExpanded.map((x, j) => (idx === j) ? !x : x)
    //     }), () => {
    //         // causes the parent to perform a (delayed) reflow
    //         this.props.onReportToggled();
    //     });
    // };

    render() {
        const variant = this.props.variant;
        const data = this.props.alleleFrequencyData;
        const exacGraph = _.find(data, function(dd) {
                                return dd.source === "ExAC";
                            }).chart[0];
        const renderedExacGraph = exacGraph.replace(variant, exacGraph.prop);
        const exacData = _.find(data, function(dd) {
                                return dd.source === "ExAC";
                            }).data;
        const renderedExacData = this.getRows("ExAC", exacData, variant);
        const thousandGenomesGraph = _.find(data, function(dd) {
                                return dd.source === "1000 Genomes";
                            }).chart[0];
        const renderedThousandGenomesGraph = thousandGenomesGraph.replace(variant, thousandGenomesGraph.prop);
        const thousandGenomesData = _.find(data, function(dd) {
                                return dd.source === "1000 Genomes";
                            }).data;
        const renderedThousandGenomesData = this.getRows("1000 Genomes", thousandGenomesData, variant);
        const espData = _.find(data, function(dd) {
                                return dd.source === "ESP";
                            }).data;
        const renderedEspData = this.getRows("ESP", espData, variant);

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
                    <Table>
                        <tbody>
                            {renderedExacGraph}
                            {renderedExacData}
                            {renderedThousandGenomesGraph}
                            {renderedThousandGenomesData}
                            {renderedEspData}
                        </tbody>
                    </Table>
                </Panel>
            </div>
        );
    };
}
