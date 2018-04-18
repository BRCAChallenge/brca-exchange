/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";
import {Panel} from 'react-bootstrap';
import {AlleleFrequencyField} from "./AlleleFrequencyField";
const _ = require('underscore');


const fieldsOfInterest = {
    'ExAC Graph': true,
    'ExAC Data': true,
    '1000 Genomes Graph': true,
    '1000 Genomes Data': true,
    'ESP Data': true
};

export default class AlleleFrequenciesTile extends React.Component {
    constructor(props) {
        super(props);

        this.state = fieldsOfInterest;

        this.fieldToggled = this.fieldToggled.bind(this);
        this.setAllFieldsExpansion = this.setAllFieldsExpansion.bind(this);
    }

    setAllFieldsExpansion(e, newExpansion) {
        e.stopPropagation();

        this.setState(function() {
            // updates all values in fieldsOfInterest to match newExpansion
            return Object.keys(_.clone(fieldsOfInterest)).forEach( function(key, idx, arr) {
                arr[key] = newExpansion;
            });
        }, () => {
            // causes the parent to perform a (delayed) reflow
            this.props.onFrequencyFieldToggled();
        });
    }

    fieldToggled(fieldName) {
        this.setState({
            [fieldName]: !this.state[fieldName]
        }, () => {
            // causes the parent to perform a (delayed) reflow
            this.props.onFrequencyFieldToggled();
        });
    }

    render() {
        const variant = this.props.variant;
        const data = this.props.alleleFrequencyData;

        const exacGraph = [_.find(data, function(dd) {
                                return dd.source === "ExAC";
                            }).chart[0], "ExAC Graph"];

        const exacData = [_.find(data, function(dd) {
                                return dd.source === "ExAC";
                            }).data, "ExAC Data"];

        const thousandGenomesGraph = [_.find(data, function(dd) {
                                return dd.source === "1000 Genomes";
                            }).chart[0], "1000 Genomes Graph"];

        const thousandGenomesData = [_.find(data, function(dd) {
                                return dd.source === "1000 Genomes";
                            }).data, "1000 Genomes Data"];

        const espData = [_.find(data, function(dd) {
                                return dd.source === "ESP";
                            }).data, "ESP Data"];

        const alleleFrequencyFields = [exacGraph, exacData, thousandGenomesGraph, thousandGenomesData, espData].map((field, idx) => {
            let fieldValue = field[0];
            let fieldName = field[1];
            let expanded = this.state[fieldName];
            return (
                <AlleleFrequencyField
                    key={idx}
                    idx={idx}
                    variant={variant}
                    field={fieldValue}
                    fieldName={fieldName}
                    expanded={expanded}
                    onFieldToggled={this.fieldToggled}
                    hideEmptyItems={this.props.hideEmptyItems}
                    showHelp={this.props.showHelp}
                />
            );
        });

        // TODO: figure out how to determine if everything is empty
        const allEmpty = false;
        // const allEmpty = _.every([renderedExacData, renderedThousandGenomesData, renderedEspData], function(data) {
        //                         return data[1] === true;
        //                     }) && !variant.Variant_in_1000_Genomes && !variant.Variant_in_ExAC;

        // create the source panel itself now
        const groupTitle = `source-panel-${this.props.sourceName}`;
        const header = (
            <h3 style={{display: 'flex', flexDirection: 'row'}}>
                <a style={{flexGrow: 1}} href="#" onClick={(event) => this.props.onChangeGroupVisibility(groupTitle, event)}>
                    {this.props.groupTitle}
                </a>

                <a title='collapse all fields'
                    onClick={(event) => this.setAllFieldsExpansion(event, false)}
                    style={{cursor: 'pointer', marginRight: '10px'}}>
                    <i className="fa fa-angle-double-up" aria-hidden="true" />
                </a>

                <a title='expand all fields'
                    onClick={(event) => this.setAllFieldsExpansion(event, true)}
                    style={{cursor: 'pointer'}}>
                    <i className="fa fa-angle-double-down" aria-hidden="true" />
                </a>
            </h3>
        );

        return (
            <div key={`group_collection-${groupTitle}`} className={ allEmpty && this.state.hideEmptyItems ? "group-empty variant-detail-group" : "variant-detail-group" }>
                <Panel
                    header={header}
                    collapsable={true}
                    defaultExpanded={localStorage.getItem("collapse-group_" + groupTitle) !== "true"}>
                    {alleleFrequencyFields}
                </Panel>
            </div>
        );
    }
};



// const AlleleFrequenciesTile = React.createClass({

//     mixins: [CollapsableMixin],

//     getInitialState: function () {
//         return {
//             // true means field is expanded, false is collapsed
//             'ExAC Graph': true,
//             'ExAC Data': true,
//             '1000 Genomes Graph': true,
//             '1000 Genomes Data': true,
//             'ESP Data': true
//         };
//     },

//     // getRowsAndDetermineIfEmpty: function(source, data, variant) {
//     //     // let styles = this.getCollapsableClassSet();
//     //     let rowsEmpty = 0;

//     //     const rows = _.map(data, (rowDescriptor, idx) => {
//     //         let {prop, title} = rowDescriptor;
//     //         let rowItem;

//     //         if (variant[prop] !== null) {
//     //             rowItem = util.getFormattedFieldByProp(prop, variant);
//     //         }

//     //         let isEmptyValue = util.isEmptyField(variant[prop]);

//     //         if (isEmptyValue) {
//     //             rowsEmpty += 1;
//     //             rowItem = '-';
//     //         }

//     //         return (
//     //             <tr key={prop} className={ (isEmptyValue && this.props.hideEmptyItems) ? "variantfield-empty" : "" }>
//     //                 { rowDescriptor.tableKey !== false &&
//     //                     <KeyInline tableKey={title} onClick={(event) => this.props.showHelp(event, title)}/>
//     //                 }
//     //                 <td colSpan={rowDescriptor.tableKey === false ? 2 : null} ><span className={"row-value" }>{rowItem}</span></td>
//     //             </tr>
//     //         );
//     //     });
//     //     const allEmpty = rowsEmpty >= data.length;

//     //     return [rows, allEmpty];
//     // },

//     getCollapsableDOMNode: function() {
//         return this.refs.panel.getDOMNode();
//     },

//     getCollapsableDimensionValue: function() {
//         return this.refs.panel.getDOMNode().scrollHeight;
//     },

//     handleToggle: function(e, fieldName) {
//         e.preventDefault();
//         this.setState({
//             [fieldName]: !this.state[fieldName]
//         }, () => {
//             // causes the parent to perform a (delayed) reflow
//             this.props.onFrequencyFieldToggled();
//         });
//     },

//     setAllFieldsExpansion: function(e, newExpansion) {
//         e.stopPropagation();

//         this.setState({
//             'ExAC Graph': newExpansion,
//             'ExAC Data': newExpansion,
//             '1000 Genomes Graph': newExpansion,
//             '1000 Genomes Data': newExpansion,
//             'ESP Data': newExpansion
//         }, () => {
//             // causes the parent to perform a (delayed) reflow
//             this.props.onFrequencyFieldToggled();
//         });
//     },

//     generateHeader: function(field, fieldName) {
//         return (
//             <div className={`allele-frequency-header ${this.state[fieldName] ? 'expanded' : ''}`} onClick={(e) => this.handleToggle(e, fieldName)}>
//                 <div className="allele-frequency-cell allele-frequency-label">
//                     {
//                         this.state[fieldName]
//                             ? <i className="fa fa-caret-down" aria-hidden="true" />
//                             : <i className="fa fa-caret-right" aria-hidden="true" />
//                     }
//                     &nbsp;
//                     <span>{fieldName}</span>
//                 </div>
//             </div>
//         );
//     },

//     generateCategoryTable: function(field, fieldName) {
//         let renderedField;
//         if (Array.isArray(field)) {
//             renderedField = field[0];
//         } else {
//             renderedField = field;
//         }

//         let styles = this.getCollapsableClassSet();

//         return (
//             <div>
//                 <div style={{marginBottom: 0, borderTop: 'solid 2px #ccc'}}>
//                 {
//                     this.generateHeader(field, fieldName)
//                 }
//                 </div>

//                 <div ref='panel' className={classNames(styles)}>
//                     <Table key={`frequency-name-${fieldName}`} id={fieldName}>
//                         <tbody>
//                             {renderedField}
//                         </tbody>
//                     </Table>
//                 </div>
//             </div>
//         );
//     },

//     render: function() {
//         const variant = this.props.variant;
//         const data = this.props.alleleFrequencyData;

//         const exacGraph = _.find(data, function(dd) {
//                                 return dd.source === "ExAC";
//                             }).chart[0];
//         const renderedExacGraph = exacGraph.replace(variant, exacGraph.prop);

//         const exacData = _.find(data, function(dd) {
//                                 return dd.source === "ExAC";
//                             }).data;
//         const renderedExacData = this.getRowsAndDetermineIfEmpty("ExAC", exacData, variant);

//         const thousandGenomesGraph = _.find(data, function(dd) {
//                                 return dd.source === "1000 Genomes";
//                             }).chart[0];
//         const renderedThousandGenomesGraph = thousandGenomesGraph.replace(variant, thousandGenomesGraph.prop);

//         const thousandGenomesData = _.find(data, function(dd) {
//                                 return dd.source === "1000 Genomes";
//                             }).data;
//         const renderedThousandGenomesData = this.getRowsAndDetermineIfEmpty("1000 Genomes", thousandGenomesData, variant);

//         const espData = _.find(data, function(dd) {
//                                 return dd.source === "ESP";
//                             }).data;
//         const renderedEspData = this.getRowsAndDetermineIfEmpty("ESP", espData, variant);

//         const allEmpty = _.every([renderedExacData, renderedThousandGenomesData, renderedEspData], function(data) {
//                                 return data[1] === true;
//                             }) && !variant.Variant_in_1000_Genomes && !variant.Variant_in_ExAC;

//         // create the source panel itself now
//         const groupTitle = `source-panel-${this.props.sourceName}`;
//         const header = (
//             <h3 style={{display: 'flex', flexDirection: 'row'}}>
//                 <a style={{flexGrow: 1}} href="#" onClick={(event) => this.props.onChangeGroupVisibility(groupTitle, event)}>
//                     {this.props.groupTitle}
//                 </a>

//                 <a title='collapse all fields'
//                     onClick={(event) => this.setAllFieldsExpansion(event, false)}
//                     style={{cursor: 'pointer', marginRight: '10px'}}>
//                     <i className="fa fa-angle-double-up" aria-hidden="true" />
//                 </a>

//                 <a title='expand all fields'
//                     onClick={(event) => this.setAllFieldsExpansion(event, true)}
//                     style={{cursor: 'pointer'}}>
//                     <i className="fa fa-angle-double-down" aria-hidden="true" />
//                 </a>
//             </h3>
//         );

//         return (
//             <div key={`group_collection-${groupTitle}`} className={ allEmpty && this.state.hideEmptyItems ? "group-empty variant-detail-group" : "variant-detail-group" }>
//                 <Panel
//                     header={header}
//                     collapsable={true}
//                     defaultExpanded={localStorage.getItem("collapse-group_" + groupTitle) !== "true"}>
//                     {this.generateCategoryTable(renderedExacGraph, 'ExAC Graph')}
//                     {this.generateCategoryTable(renderedExacData, 'ExAC Data')}
//                     {this.generateCategoryTable(renderedThousandGenomesGraph, '1000 Genomes Graph')}
//                     {this.generateCategoryTable(renderedThousandGenomesData, '1000 Genomes Data')}
//                     {this.generateCategoryTable(renderedEspData, 'ESP Data')}
//                 </Panel>
//             </div>
//         );
//     }
// });

// module.exports = AlleleFrequenciesTile;
