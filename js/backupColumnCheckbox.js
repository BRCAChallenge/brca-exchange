/*global module: false, require: false */
'use strict';

var React = require('react');
var PureRenderMixin = require('./PureRenderMixin'); // deep-equals version of PRM
var {Input} = require('react-bootstrap');

var ColumnCheckbox = React.createClass({
    mixins: [PureRenderMixin],
    onChange: function (e) {
        return this.props.onChange(e.target.value);
    },
    render: function () {
        var {label, title, initialCheck} = this.props;
        var subColumnLabel = '';
        if (title === 'Gene') {
            subColumnLabel = 'ENIGMA Columns:';
        } else if (title === 'Allele frequency (1000 Genomes)') {
            subColumnLabel = '1000 Genomes Columns:';
        } else if (title === 'Allele origin (ClinVar)') {
            subColumnLabel = 'ClinVar Columns:';
        }
        return (
            <div>
            <div><h5>{subColumnLabel}</h5></div>
            <Input type="checkbox" label={title} checked={initialCheck[label].selectVal} onChange={this.onChange} />
            </div>
        );
    }
});

module.exports = ColumnCheckbox;
