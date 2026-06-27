/*eslint-env browser */
'use strict';

import React from 'react';
import slugify from './slugify';
import { withRouter } from 'react-router-dom';
import AutoSuggest from 'react-autosuggest';
import _ from 'underscore';
import $ from 'jquery';
import config from './config';

import './css/Autosuggest.css';

function getSuggestions(value, callback, release) {
    const matchStr = encodeURIComponent(value.toLowerCase());
    let suggestionsEndpoint = `${config.backend_url}/data/suggestions/?term=${matchStr}`;
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

class VariantSearchInner extends React.Component {
    static defaultProps = { onSearch: () => {} };

    constructor(props) {
        super(props);
        this.state = {
            value: props.value || '',
            release: props.release,
            placeholder: 'search for "c.1105G>A", "brca1" or "IVS7+1037T>C"',
            suggestions: []
        };

        this.onClick = this.onClick.bind(this);
        this.onClickSearchButton = this.onClickSearchButton.bind(this);
        this.showHelp = this.showHelp.bind(this);
        this.onChange = this.onChange.bind(this);
        this.onSubmit = this.onSubmit.bind(this);
        this.onFocus = this.onFocus.bind(this);
        this.onBlur = this.onBlur.bind(this);
        this.onFetchSuggestions = this.onFetchSuggestions.bind(this);
        this.onClearSuggestions = this.onClearSuggestions.bind(this);
    }

    onClick(ev) {
        ev.stopPropagation();
        this.props.onSearch(this.state.value);
    }

    onClickSearchButton() {
        this.props.onSearch(this.state.value);
    }

    showHelp(title) {
        this.props.history.push(`/help#${slugify(title)}`);
    }

    onChange(event, { newValue: value }) {
        const { onChange } = this.props;
        // Avoid loops when props update triggers onChange
        if (value !== this.state.value && onChange) {
            onChange(value);
        }
        this.setState({ value: value || '', release: this.props.release });
    }

    onSubmit(ev) {
        ev.preventDefault();
        this.props.onSearch(this.state.value);
    }

    onFocus() {
        this.setState({ placeholder: '' });
    }

    onBlur() {
        this.setState({ placeholder: 'search for "c.1105G>A", "brca1" or "IVS7+1037T>C"' });
    }

    onFetchSuggestions({ value }) {
        getSuggestions(
            value,
            (_error, results) => {
                this.setState({ suggestions: results });
            },
            this.state.release
        );
    }

    onClearSuggestions() {
        this.setState({ suggestions: [] });
    }

    componentDidUpdate(prevProps) {
        if (prevProps.value !== this.props.value) {
            this.setState({ value: this.props.value || '' });
        }
    }

    render() {
        const { id, onSearch } = this.props;
        const { value, suggestions } = this.state;

        return (
            <div className='search-box'>
                <form onSubmit={this.onSubmit} style={{ display: 'inline' }}>
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
                        />

                        <span
                            className="fa fa-search search-box-icon"
			    role="button"
			    tabIndex={0}
			    aria-label="Search"
                            onClick={this.onClickSearchButton}
                        />
                    </div>
                </form>
            </div>
        );
    }
}

export default withRouter(VariantSearchInner);

