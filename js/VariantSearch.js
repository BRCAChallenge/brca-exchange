/*global require: false, module: false */
'use strict';

var React = require('react');
var slugify = require('./slugify');
var {Button} = require('react-bootstrap');
var {Navigation} = require('react-router');

var VariantSearch = React.createClass({
	mixins: [Navigation],
	onClick: function (ev) {
		ev.stopPropagation();
		var value = React.findDOMNode(this.refs.input).value;
		this.props.onSearch(value);
	},
	onKeyDown: function (ev) {
		if (ev.key === 'Enter') {
			this.props.onSearch(ev.target.value);
		}
	},
	showHelp: function (title) {
		this.transitionTo(`/help#${slugify(title)}`);
	},
	onChange: function (ev) {
		var {onChange} = this.props;
		if (onChange) {
			onChange(ev.target.value);
		}
	},
	render: function () {
		var {value} = this.props;
		return (
			<div className='search-box help-target'>
				<input ref='input' value={value} onChange={this.onChange} onKeyDown={this.onKeyDown} placeholder="Search Variant such as c.1105G>A"></input>
				<span className='text-nowrap'>
					<Button onClick={this.onClick} className='btn-xs'>
						<span>&nbsp;&nbsp;</span>
						<span className="glyphicon glyphicon-search"></span>
					</Button>
					<span onClick={() => this.showHelp('Searching')}
						className="glyphicon glyphicon-question-sign superscript help"/>
				</span>
			</div>
		);
	}
});

module.exports = VariantSearch;
