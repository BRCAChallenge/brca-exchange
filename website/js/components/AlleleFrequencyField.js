/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";
import {CollapsableMixin, Table} from "react-bootstrap";
import classNames from 'classnames';
import util from '../util';
import KeyInline from './KeyInline';
import slugify from "../slugify";
const _ = require('underscore');


const AlleleFrequencyField = React.createClass({
    mixins: [CollapsableMixin],

    getCollapsableDOMNode: function() {
        return this.refs.panel.getDOMNode();
    },

    getCollapsableDimensionValue: function() {
        return this.refs.panel.getDOMNode().scrollHeight;
    },

    handleToggle: function(e, fieldName) {
        e.preventDefault();

        // ask our parent to toggle us
        this.props.onFieldToggled(fieldName);
    },

    generateHeader: function(field, fieldName) {
        return (
            <div className={`allele-frequency-header ${this.props.expanded ? 'expanded' : ''}`} onClick={(e) => this.handleToggle(e, fieldName)}>
                <div className="allele-frequency-cell allele-frequency-label">
                    {
                        this.props.expanded
                            ? <i className="fa fa-caret-down" aria-hidden="true" />
                            : <i className="fa fa-caret-right" aria-hidden="true" />
                    }
                    &nbsp;
                    <span>{fieldName}</span>
                </div>
            </div>
        );
    },

    getRowsAndDetermineIfEmpty(source, data, variant) {
        let rowsEmpty = 0;

        const rows = _.map(data, (rowDescriptor) => {
            let {prop, title} = rowDescriptor;
            let rowItem;

            if (variant[prop] !== null) {
                rowItem = util.getFormattedFieldByProp(prop, variant);
            }

            let isEmptyValue = util.isEmptyField(variant[prop]);

            if (isEmptyValue) {
                rowsEmpty += 1;
                rowItem = '-';
            }

            return (
                <tr key={prop} className={ (isEmptyValue && this.props.hideEmptyItems) ? "variantfield-empty" : "" }>
                    { rowDescriptor.tableKey !== false &&
                        (
                            <KeyInline tableKey={title}
                                tooltip={this.props.tooltips && this.props.tooltips[slugify(prop)]}
                                onClick={(event) => this.props.showHelp(event, prop)}
                            />
                        )
                    }
                    <td colSpan={rowDescriptor.tableKey === false ? 2 : null} ><span className={"row-value" }>{rowItem}</span></td>
                </tr>
            );
        });
        const allEmpty = rowsEmpty >= data.length;

        return [rows, allEmpty];
    },

    render: function() {
        const {field, fieldName, variant, hideEmptyItems} = this.props;
        let renderedRows;
        let allEmpty = false;
        let styles = this.getCollapsableClassSet();
        let isChart = false;

        if (fieldName === "ExAC (Graphical)") {
            renderedRows = field.replace(variant, field.prop);
            if (!variant.Variant_in_ExAC) {
                allEmpty = true;
            }
            isChart = true;
        } else if (fieldName === "ExAC (Numerical)") {
            renderedRows = this.getRowsAndDetermineIfEmpty("ExAC", field, variant);
        } else if (fieldName === "1000 Genomes (Graphical)") {
            renderedRows = field.replace(variant, field.prop);
            if (!variant.Variant_in_1000_Genomes) {
                allEmpty = true;
            }
            isChart = true;
        } else if (fieldName === "1000 Genomes (Numerical)") {
            renderedRows = this.getRowsAndDetermineIfEmpty("1000 Genomes", field, variant);
        } else if (fieldName === "ESP (Numerical)") {
            renderedRows = this.getRowsAndDetermineIfEmpty("ESP", field, variant);
        }

        if (Array.isArray(renderedRows)) {
            allEmpty = renderedRows[1];
            renderedRows = renderedRows[0];
        }

        return (
            <div className={ allEmpty && (isChart || hideEmptyItems) ? "group-empty" : "" }>
                <div style={{marginBottom: 0, borderTop: 'solid 2px #ccc'}}>
                {
                    this.generateHeader(field, fieldName)
                }
                </div>

                <div ref='panel' className={allEmpty && isChart ? "group-empty" : classNames(styles)}>
                    <Table key={`allele-frequency-name-${fieldName}`}>
                        <tbody>
                        { renderedRows }
                        </tbody>
                    </Table>
                </div>
            </div>
        );
    }
});

module.exports = {
    AlleleFrequencyField,
};
