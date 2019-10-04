'use strict';

import React from 'react';
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
const SearchController = React.createClass({
    getInitialState() {
        return {
            searchTerm: '',
            searching: false,
            matched: 0,
            currentMark: null
        };
    },

    searchChanged(e) {
        this.setState({
            searchTerm: e.target.value,
            searching: true
        }, () => {
            this.debouncedSearchResponse();
        });
    },

    searchResponse() {
        if (!this.state.searchTerm || this.state.searchTerm === '') {
            // remove any marks if they cleared the search
            this.setState({
                matched: 0,
                currentMark: null,
                searching: false
            }, () => {
                this.searcher.unmark();

                $(this.props.target).find('*[data-expander-id]').each((idx, elem) => {
                    this.props.setExpansion($(elem).data('expander-id'), false);
                });
            });
            return;
        }

        // perform full matching against 'target'
        // (first we unmark, then mark, then deal with the match results)
        this.searcher.unmark({
            done: () => {
                this.pendingUpdate = new Set();

                // then iteratively expand while searching for marks
                this.searcher.mark(this.state.searchTerm, {
                    element: 'span',
                    className: 'highlighted',
                    done: (totalMarks) => {
                        this.setState({
                            searching: false,
                            currentMark: null,
                            matched: totalMarks
                        });

                        // set all the elements that can be toggled to their match status
                        $('*[data-expander-id]').each((idx, elem) => {
                            const targetID = $(elem).data('expander-id');
                            this.props.setExpansion(targetID, this.pendingUpdate.has(targetID));
                        });
                    },
                    each: (elem) => {
                        const me = this;

                        // check if it has ancestors that need to be expanded and add them to the expanded list
                        $(elem)
                            .click(function() {
                                const $highlightSet = $('.highlighted').removeClass("focused");

                                // set this element to the currently-selected index
                                $(this).addClass("focused");
                                me.setState({ currentMark: $highlightSet.index(this) });
                            })
                            .parents('*[data-expander-id]').each((idx, elem) => {
                                this.pendingUpdate.add($(elem).data('expander-id'));
                            });
                    }
                });
            }
        });
    },

    componentWillReceiveProps(nextProps) {
        if (nextProps.researchMode !== this.props.researchMode) {
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
    },

    componentWillMount() {
        this.searcher = new Mark(this.props.target);
        this.debouncedSearchResponse = debounce(this.searchResponse, 300);

        this.navForward = this.navMarks.bind(this, true);
        this.navBackward = this.navMarks.bind(this, false);
    },

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
    },

    render() {
        return (
            <div className="input-group has-feedback has-search">
                <div className="input-group-addon">
                    <span className={`glyphicon ${this.state.searching ? "glyphicon-refresh glyphicon-spin" : "glyphicon-search"}`} />
                </div>
                <input type="text" className="form-control" placeholder="Search" value={this.state.searchTerm} onChange={this.searchChanged} />
                {
                    (this.state.matched > 0) && (
                        <span className="input-group-addon">
                         { this.state.currentMark !== null && `${this.state.currentMark + 1} / ` }
                            { this.state.matched}
                        </span>
                    )
                }
                <div className="input-group-btn">
                    <button type="button" disabled={this.state.matched <= 0} onClick={this.navForward} className="btn btn-default">
                        <span className="glyphicon glyphicon-triangle-bottom" />
                    </button>
                    <button type="button" disabled={this.state.matched <= 0} onClick={this.navBackward} className="btn btn-default">
                        <span className="glyphicon glyphicon-triangle-top" />
                    </button>
                </div>
            </div>
        );
    }
});
SearchController.propTypes = {
    target: React.PropTypes.string.isRequired,
    setExpansion: React.PropTypes.func.isRequired
};

export default SearchController;

