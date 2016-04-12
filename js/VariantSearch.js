/*global require: false, module: false */
/*eslint-env browser */
'use strict';

var React = require('react');
var slugify = require('./slugify');
var {Navigation} = require('react-router');
var AutoSuggest = require('react-autosuggest');
var _ = require('underscore');
var $ = require('jquery');
var config = require('./config')

require('./css/Autosuggest.css');


function getSuggestions(value, callback) {
    var matchStr = encodeURIComponent(value.toLowerCase());
    var suggestionsEndpoint = `${config.backend_url}/data/suggestions/?term=${matchStr}`;
    $.ajax({
        url: suggestionsEndpoint,
        dataType: 'json',
        success: function (data) {
            var suggestions = _.flatten(_.values(data.suggestions));
            setTimeout(() => callback(null, suggestions), 300);
        }.bind(this),
        error: function (xhr, status, err) {
            callback(new Error("Couldn't get suggestions"));
        }.bind(this)
    });
}

function renderSuggestion(suggestion, input) {
    var maxLengthToDisplay = 50;
    return (
        <span>
            <strong>{suggestion.slice(0, input.length)}</strong>
            {suggestion.slice(input.length, maxLengthToDisplay)}
            {(suggestion.length > maxLengthToDisplay) ? "..." : ""}
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
        // XXX We're getting an onChange event when props are updated, which
        // leads to a loop of sorts. Check if value has actually changed before
        // calling.
        if (value !== this.state.value && onChange) {
            onChange(value);
        }
        this.setState({value: value});
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
                    <input type='submit' className='input-sm'style={{display: 'none'}} />
                    <div className='text-nowrap help-target'>
                        <AutoSuggest
                            id={id}
                            className='dropdown open'
                            cache={false}
                            value={value}
                            inputAttributes={{
                                className: 'variant-search-input',
                                placeholder: "search for \"c.1105G>A\", \"brca1\" or \"IVS7+1037T>C\"",
                                onChange: this.onChange
                            }}
                            showWhen={input => input.trim().length > 0}
                            suggestions={_.debounce(getSuggestions, 300, true)}
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
