/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";
import {CollapsableMixin, Table} from "react-bootstrap";
import classNames from 'classnames';
import util from '../util';
import KeyInline from './KeyInline';
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
                        <KeyInline tableKey={title} onClick={(event) => this.props.showHelp(event, title)}/>
                    }
                    <td colSpan={rowDescriptor.tableKey === false ? 2 : null} ><span className={"row-value" }>{rowItem}</span></td>
                </tr>
            );
        });
        const allEmpty = rowsEmpty >= data.length;

        return [rows, allEmpty];
    },

    render: function() {
        const {field, fieldName, variant} = this.props;
        let renderedRows;

        if (fieldName === "ExAC Graph") {
            renderedRows = field.replace(variant, field.prop);
        } else if (fieldName === "ExAC Data") {
            renderedRows = this.getRowsAndDetermineIfEmpty("ExAC", field, variant);
        } else if (fieldName === "1000 Genomes Graph") {
            renderedRows = field.replace(variant, field.prop);
        } else if (fieldName === "1000 Genomes Data") {
            renderedRows = this.getRowsAndDetermineIfEmpty("1000 Genomes", field, variant);
        } else if (fieldName === "ESP Data") {
            renderedRows = this.getRowsAndDetermineIfEmpty("ESP", field, variant);
        }

        let styles = this.getCollapsableClassSet();

        return (
            <div>
                <div style={{marginBottom: 0, borderTop: 'solid 2px #ccc'}}>
                {
                    this.generateHeader(field, fieldName)
                }
                </div>

                <div ref='panel' className={classNames(styles)}>
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
