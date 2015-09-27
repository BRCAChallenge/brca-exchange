/*global require: false, module: false */
/*eslint-env browser */
'use strict';

var React = require('react');
var slugify = require('./slugify');
var {Navigation} = require('react-router');
var AutoSuggest = require('react-autosuggest');
var _ = require('underscore');
require('./css/Autosuggest.css');

var maxSuggestions = 10;
function getSuggestions(data, input) {
	var matchStr = input.toLowerCase();
	return _.first(_.filter(data, s => s.toLowerCase().indexOf(matchStr) === 0),
			maxSuggestions);
}

function renderSuggestion(suggestion, input) {
	return (
		<span>
			<strong>{suggestion.slice(0, input.length)}</strong>
			{suggestion.slice(input.length)}
		</span>
   );
}

var VariantSearch = React.createClass({
	mixins: [Navigation],
	onClick: function (ev) {
		ev.stopPropagation();
		var value = React.findDOMNode(this.refs.input).value;
		this.props.onSearch(value);
	},
	showHelp: function (title) {
		this.transitionTo(`/help#${slugify(title)}`);
	},
	onChange: function (value) {
		var {onChange} = this.props;
		if (onChange) {
			onChange(value);
		}
		this.setState({value: value});
	},
	suggest: function (input, callback) {
		var {suggestions} = this.props;
		// Invoke asynchronously. This makes more sense if doing an ajax call.
		this.cb = setTimeout(() =>
				callback(null, getSuggestions(suggestions, input, callback), 0));
	},
	componentWillUnmount: function () {
		clearTimeout(this.cb);
	},
	onSubmit: function (ev) {
		ev.preventDefault();
		this.props.onSearch(this.state.value);
	},
	getDefaultProps: function () {
		return {
			onSearch: () => {}
		};
	},
	getInitialState: function () {
		return {
			value: this.props.value
		};
	},
	componentWillReceiveProps: function (newProps) {
		this.setState({value: newProps.value});
	},
	render: function () {
		var {id, onSearch} = this.props,
			{value} = this.state;

		return (
			<div className='search-box'>
				<form onSubmit={this.onSubmit} style={{display: 'inline'}}>
					<input type='submit' style={{display: 'none'}} />
					<div className='text-nowrap help-target'>
						<AutoSuggest
							id={id}
							className='dropdown open'
							cache={false}
							value={value}
							inputAttributes={{
								className: 'variant-search-input',
								placeholder: "Search Variant such as c.1105G>A",
								onChange: this.onChange
							}}
							showWhen={input => input.trim().length > 1}
							suggestions={this.suggest}
							onSuggestionSelected={v => onSearch(v)}
							suggestionRenderer={renderSuggestion}
							ref='input' />
						<span onClick={() => this.showHelp('Searching')}
							className="glyphicon glyphicon-question-sign superscript help"/>
					</div>
				</form>
			</div>
		);
	}
});

module.exports = VariantSearch;
