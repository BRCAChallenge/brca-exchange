/*global module: false, require: false */
'use strict';

var React = require('react');
var PureRenderMixin = require('./PureRenderMixin'); // deep-equals version of PRM
var {Input} = require('react-bootstrap');
var _ = require('underscore');

var SelectField = React.createClass({
    mixins: [PureRenderMixin],
    onChange: function (e) {
        return this.props.onChange(e.target.value);
    },
    render: function () {
        var {options, label, value} = this.props,
            opels = _.map(options, v => <option key={v} value={v}>{v}</option>);

        return (
            <Input type="select" value={value} label={label} onChange={this.onChange}>
                {opels}
            </Input>
        );
    }
});

module.exports = SelectField;
