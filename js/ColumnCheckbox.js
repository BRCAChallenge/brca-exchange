/*global module: false, require: false */
'use strict';

var React = require('react');
var PureRenderMixin = require('./PureRenderMixin'); // deep-equals version of PRM
var {Input} = require('react-bootstrap');
var _ = require('underscore');

var ColumnCheckbox = React.createClass({
	mixins: [PureRenderMixin],
	onChange: function (e) {
		return this.props.onChange(e.target.value);
	},
	render: function () {
		var {options, label, prop} = this.props,
			opels = _.map(options, v => <option key={v} value={v}>{v}</option>);
		return (
            <div>
			<Input type="checkbox" label={label} onChange={this.onChange}>
			</Input>
            </div>
		);
	}
});

module.exports =ColumnCheckbox;
