/*global require: false, module: false */
/*eslint-env browser */
'use strict';

import React from 'react';
import ReactDOM from 'react-dom';
import slugify from './slugify';
import {Navigation} from 'react-router';
import AutoSuggest from 'react-autosuggest';
import * as _ from 'underscore';
import * as $ from 'jquery';
import config from './config';

import './css/Autosuggest.css';


function getSuggestions(value, callback, release) {
    var matchStr = encodeURIComponent(value.toLowerCase());
    var suggestionsEndpoint = `${config.backend_url}/data/suggestions/?term=${matchStr}`;
    // If a release is specified, include it in the request
    if (release) {
        suggestionsEndpoint += `&release=${release}`;
    }
    $.ajax({
        url: suggestionsEndpoint,
        dataType: 'json',
        success: function (data) {
            callback(null, _.flatten(_.values(data.suggestions)));
        },
        error: function () {
            callback(new Error("Couldn't get suggestions"));
        }
    });
}

function renderSuggestion(suggestion, { query }) {
    const maxLengthToDisplay = 50;
    return (
        <span>
            <strong>{suggestion.slice(0, query.length)}</strong>
            {suggestion.slice(query.length, maxLengthToDisplay)}
            {(suggestion.length > maxLengthToDisplay) ? "..." : ""}
        </span>
    );
}

var VariantSearch = React.createClass({
    mixins: [Navigation],
    onClick: function (ev) {
        ev.stopPropagation();
        var value = ReactDOM.findDOMNode(this.refs.input).value;
        this.props.onSearch(value);
    },

    onClickSearchButton: function () {
        this.props.onSearch(this.state.value);
    },

    showHelp: function (title) {
        this.transitionTo(`/help#${slugify(title)}`);
    },

    onChange: function (event, { newValue: value }) {
        var {onChange} = this.props;
        // XXX We're getting an onChange event when props are updated, which
        // leads to a loop of sorts. Check if value has actually changed before
        // calling.
        if (value !== this.state.value && onChange) {
            onChange(value);
        }
        this.setState({value: value, release: this.props.release});
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
            value: this.props.value || '',
            release: this.props.release,
            placeholder: "search for \"c.1105G>A\", \"brca1\" or \"IVS7+1037T>C\"",
            suggestions: this.props.suggestions || []
        };
    },

    onFocus: function () {
        this.setState({placeholder: ""});
    },

    onBlur: function() {
        this.setState({placeholder: "search for \"c.1105G>A\", \"brca1\" or \"IVS7+1037T>C\""});
    },

    onFetchSuggestions({ value }) {
        getSuggestions(value, (error, results) => {
            this.setState({ suggestions: results });
        }, this.state.release);
    },

    onClearSuggestions() {
        this.setState({ suggestions: [] });
    },

    componentWillReceiveProps: function (newProps) {
        this.setState({value: newProps.value});
    },

    render: function () {
        const {id, onSearch} = this.props;
        const {value, suggestions} = this.state;

        return (
            <div className='search-box'>
                <form onSubmit={this.onSubmit} style={{display: 'inline'}}>
                    <div className='text-nowrap help-target'>
                        <AutoSuggest
                            id={id}
                            className='dropdown open'
                            inputProps={{
                                className: 'variant-search-input',
                                placeholder: this.state.placeholder,
                                onChange: this.onChange,
                                onFocus: this.onFocus,
                                onBlur: this.onBlur,
                                value: value
                            }}
                            shouldRenderSuggestions={input => input.trim().length > 0}
                            onSuggestionsFetchRequested={_.debounce(this.onFetchSuggestions, 200)}
                            onSuggestionsClearRequested={this.onClearSuggestions}
                            getSuggestionValue={(x) => x}
                            suggestions={suggestions}
                            onSuggestionSelected={(event, { suggestionValue }) => onSearch(suggestionValue)}
                            renderSuggestion={renderSuggestion}
                            ref='input'
                        />

                        <span
                            className="glyphicon glyphicon-search search-box-icon"
                            onClick={this.onClickSearchButton}
                        />
                    </div>
                </form>
            </div>
        );
    },
});

module.exports = VariantSearch;
