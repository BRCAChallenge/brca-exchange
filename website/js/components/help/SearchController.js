'use strict';

import React from 'react';
import PropTypes from 'prop-types';
import Mark from "mark.js";
import {debounce} from "lodash";
import {navbarHeight} from "../../Help";
const $ = require('jquery');

export const EXTRA_SEARCH_PADDING = 8;

/**
 * The SearchController renders as a sticky search box with extra controls for navigating the search results.
 *
 * When beginning a search, all collapsible elements in the targeted element's tree will be collapsed. Strings that
 * match the query will be wrapped in <span class="highlighted"> tags, and all ancestors of the string that are
 * collapsible will be expanded.
 *
 * When the site's mode is toggled, search results will be cleared since they won't be valid for the new text anyway.
 */
class SearchController extends React.Component {
    constructor(props) {
        super(props);
        this.state = { searchTerm: '', searching: false, matched: 0, currentMark: null };
        this.searcher = new Mark(this.props.target);
        this.debouncedSearchResponse = debounce(this.searchResponse.bind(this), 300);
        this.searchChanged = this.searchChanged.bind(this);
        this.navMarks = this.navMarks.bind(this);
        this.navForward = this.navMarks.bind(this, true);
        this.navBackward = this.navMarks.bind(this, false);
    }

    searchChanged(e) {
        this.setState({
            searchTerm: e.target.value,
            searching: true
        }, () => {
            this.debouncedSearchResponse();
        });
    }

searchResponse() {
    if (!this.state.searchTerm || this.state.searchTerm === '') {
        this.setState({ matched: 0, currentMark: null, searching: false });
        this.searcher.unmark();
        $('*[data-expander-id]').each((idx, elem) => {
            this.props.setExpansion($(elem).data('expander-id'), false);
        });
        return;
    }

    // PASS 1: find which cards contain matches, expand them, but don't keep marks
    this.searcher.unmark({
        done: () => {
            const expandSet = new Set();

            this.searcher.mark(this.state.searchTerm, {
                element: 'span',
                className: 'highlighted',
                each: (elem) => {
                    $(elem).parents('*[data-expander-id]').each((idx, parent) => {
                        expandSet.add($(parent).data('expander-id'));
                    });
                },
                done: () => {
                    // expand cards with matches, collapse those without
                    $('*[data-expander-id]').each((idx, elem) => {
                        const id = $(elem).data('expander-id');
                        this.props.setExpansion(id, expandSet.has(id));
                    });

                    // PASS 2: after React re-renders and Collapse animates, re-apply marks
                    setTimeout(() => {
                        this.searcher.unmark({
                            done: () => {
                                this.searcher.mark(this.state.searchTerm, {
                                    element: 'span',
                                    className: 'highlighted',
                                    each: (elem) => {
                                        $(elem).click(function() {
                                            const $highlightSet = $('.highlighted').removeClass("focused");
                                            $(this).addClass("focused");
                                            this.setState({ currentMark: $highlightSet.index(this) });
                                        }.bind(this));
                                    },
                                    done: (totalMarks) => {
                                        this.setState({
                                            searching: false,
                                            currentMark: null,
                                            matched: totalMarks
                                        });
                                    }
                                });
                            }
                        });
                    }, 400); // wait for Collapse animation + React re-render
                }
            });
        }
    });
}

    componentDidUpdate(prevProps) {
        if (prevProps.researchMode !== this.props.researchMode) {
            // reinitialize if they've switched help page modes
            this.setState({
                searchTerm: '',
                searching: false,
                matched: 0,
                currentMark: null
            });

            // clear any marks and recreate the marker
            this.searcher.unmark({
                done: () => {
                    this.searcher = new Mark(this.props.target);
                }
            });
        }
    }

    navMarks(forward) {
        const $highlightSet = $('.highlighted').removeClass("focused");
        let nextMark = this.state.currentMark;

        if (nextMark === null) {
            // initialize currentMark if we haven't navigated anything previously
            nextMark = forward ? 0 : $highlightSet.length - 1;
        } else {
            // apply navigation
            nextMark = (nextMark + (forward ? 1 : -1)) % $highlightSet.length;
            if (nextMark < 0) {
                nextMark = $highlightSet.length + nextMark;
            }
        }

        this.setState({
            currentMark: nextMark
        }, () => {
            const $targetElem = $($highlightSet.get(this.state.currentMark)).addClass("focused");
            // move to whatever we navigated to
            window.scrollTo({
                // we want the element to not be covered by the navbar or the sticky search header, with some
                // extra cosmetic padding, EXTRA_PADDING, past the header as well
                top: $targetElem.offset().top - navbarHeight - (this.props.headerElem.outerHeight() + EXTRA_SEARCH_PADDING),
                behavior: 'smooth'
            });
        });
    }

    render() {
        return (
            <div className="input-group has-feedback has-search">
		<span className="input-group-text">
                    <span className={`fa ${this.state.searching ? "fa-refresh fa-spin" : "fa-search"}`} />
                </span>
                <input type="text" className="form-control" placeholder="Search" value={this.state.searchTerm} onChange={this.searchChanged} />
                {
                    (this.state.matched > 0) && (
			<span className="input-group-text"> 
                            { this.state.currentMark !== null && `${this.state.currentMark + 1} / ` }
                            { this.state.matched}
			</span>
                    )
                }
                <button type="button" disabled={this.state.matched <= 0} onClick={this.navForward} className="btn btn-default">
                    <span className="fa fa-caret-down" />
                </button>
                <button type="button" disabled={this.state.matched <= 0} onClick={this.navBackward} className="btn btn-default">
                    <span className="fa fa-caret-up" />
                </button>
            </div>
        );
    }
}
SearchController.propTypes = {
    target: PropTypes.string.isRequired,
    setExpansion: PropTypes.func.isRequired,
    researchMode: PropTypes.any,
    headerElem: PropTypes.object
};

export default SearchController;

