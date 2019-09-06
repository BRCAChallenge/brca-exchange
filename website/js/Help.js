/*eslint-env browser */
'use strict';

import classNames from "classnames";
import update from 'immutability-helper';
import BetaTag from "./components/BetaTag";
import {debounce} from "lodash";
import Mark from 'mark.js';
// import debounce from 'lodash/debounce';
const React = require('react');
const RawHTML = require('./RawHTML');
const {Grid, Col, Row, Panel, ListGroup, ListGroupItem, Glyphicon, CollapsableMixin, BootstrapMixin} = require('react-bootstrap');
const {State} = require('react-router');

var $ = require('jquery');
const slugify = require('./slugify');
const content = require('./content');

const navbarHeight = 70;
// extra padding when scrolling to a search result
const EXTRA_SEARCH_PADDING = 8;

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
                        // check if it has ancestors that need to be expanded and add them to the expanded list
                        $(elem).parents('*[data-expander-id]').each((idx, elem) => {
                            this.pendingUpdate.add($(elem).data('expander-id'));
                        });
                    }
                });
            }
        });
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

const CollapsableListItem = React.createClass({
    mixins: [State, BootstrapMixin, CollapsableMixin],

    onClick(e) {
        this.props.setExpansion(this.props.id);
        // this.setState({ expanded: !this.state.expanded });
        e.preventDefault();
    },

    getCollapsableDimensionValue() {
        return React.findDOMNode(this.refs.content).scrollHeight;
    },

    getCollapsableDOMNode() {
        if (!this.isMounted() || !this.refs || !this.refs.content) {
            return null;
        }

        return React.findDOMNode(this.refs.content);
    },

    render: function() {
        let {header, id, ...rest} = this.props;
        header = (
            <h4>
                <a href="#" onClick={this.onClick}
                    style={{color: "inherit", textDecoration: "none"}}
                    id={id} >
                    <small><Glyphicon  glyph={this.state.expanded ? "chevron-down" : "chevron-right"} /> </small>
                    <span style={{verticalAlign: "text-bottom"}}>{header}</span>
                </a>
            </h4>);
        return (
            <ListGroupItem header={header} {...rest}>
                <div className={classNames(this.getCollapsableClassSet("collapse"))}
                    ref="content">
                    { this.props.children }
                </div>
            </ListGroupItem>
        );
    }
});

const Help = React.createClass({
    mixins: [State],

    getInitialState() {
        let help = localStorage.getItem("research-mode") === 'true' ? content.helpContentResearch : content.helpContentDefault;
        const {fragment, fragRegex} = this.fragmentMatchers();

        const collapsedItems = {};

        help.forEach(({tiles}) => {
            tiles.forEach(({name, id, contents, list}) => {
                collapsedItems[id ? id : slugify(name)] = this.shouldBeExpanded(fragment, fragRegex, {name, id, contents, list});
                if (list) {
                    list.forEach(({name, id, contents}) => {
                        collapsedItems[id ? id : slugify(name)] = this.shouldBeExpanded(fragment, fragRegex, {name, id, contents});
                    });
                }
            });
        });

        return {
            searchTerm: '',
            collapsedItems
        };
    },

    fragmentMatchers() {
        const fragment =  slugify(window.location.hash.slice(1));
        return {
            fragment,
            // looks for the fragment within id attributes wrapped in quotes (in dev) or without quotes (minified)
            // we compile it once here and pass the same regex to each call to shouldBeExpanded() to save a little time
            fragRegex: new RegExp(`(id=${fragment}|id="${fragment}")`)
        };
    },

    componentWillMount() {
        // debouncing committedSearch prevents the collapsing panels from reacting too quickly to user input
        // (e.g., opening a panel immediately on 'b', then closing it on 'br', etc.)
        this.debouncedCommitSearch = debounce(() => {
            this.setState((pstate) => ({
                committedSearch: pstate.searchTerm
            }));
        }, 300);
    },

    componentDidMount() {
        const fragment = slugify(window.location.hash.slice(1));
        if (fragment !== '') {
            setTimeout(function () {
                var el = document.getElementById(fragment);
                if (el) {
                    window.scrollTo(0, el.getBoundingClientRect().top - navbarHeight);
                }
            }, 0);
        }
    },

    shouldBeExpanded(fragment, fragRegex, {name, id, contents, list}) {
        // if fragment is empty, fragRegex will be (id=|id=""), which will match anything containing 'id='
        // if there isn't a fragment, then we should never expand anything by default anyway, so just bail w/false
        if (!fragment) {
            return !!list;
        }

        let slug = id || slugify(name);
        if (slug === fragment) {
            return true;
        }
        if (contents && fragRegex.test(contents)) {
            return true;
        }
        if (list && list.some(elem => this.shouldBeExpanded(fragment, fragRegex, elem))) {
            return true;
        }
        return false;
    },

    searchChanged(event) {
        this.setState({
            searchTerm: event.target.value
        }, () => {
            this.debouncedCommitSearch();
        });
    },

    setExpansion(id, forced = null) {
        this.setState((pstate) => ({
            collapsedItems: update(pstate.collapsedItems, {
                [id]: (forced !== null) ? {$set: forced} : (visible) => visible ? !visible : true
            })
        }));
    },

    render() {
        let help = localStorage.getItem("research-mode") === 'true' ? content.helpContentResearch : content.helpContentDefault;

        const {fragment} = this.fragmentMatchers();

        var helpTiles = help.map(({section, tiles}) =>
            [<h1>{section}</h1>, tiles.map(({name, id, contents, list, reference, isBeta}) => {
                let header = [<span id={id || slugify(name)}>{name}</span>];
                // if the user clicks a reference link in a tile header, don't toggle the tile, and open the link.
                const _this = this;
                let onSelect = function (e) {
                    _this.setExpansion(e.target.getAttribute('id'));

                    if (e.target.classList.contains("help-reference-link")) {
                        e.selected = false;
                    } else {
                        e.preventDefault();
                    }
                };

                const actualId = id ? id : slugify(name);
                let body = [];
                if (contents) {
                    body.push(<RawHTML html={contents} />);
                }

                if (list) {
                    body.push(
                        <ListGroup fill>
                            {
                                list.map(({name, id, contents}) => {
                                    const localId = id ? id : slugify(name);
                                    return (
                                        <CollapsableListItem
                                            id={localId} setExpansion={this.setExpansion}
                                            expanded={this.state.collapsedItems[localId]}
                                            data-expander-id={localId}
                                            header={name}>
                                            <RawHTML html={contents} />
                                        </CollapsableListItem>
                                    );
                                })
                            }
                        </ListGroup>
                    );
                }

                if (reference) {
                    header.push(
                        <small>&nbsp;
                            <a href={reference} target="_blank">
                                <Glyphicon glyph="link" className="help-reference-link" />
                            </a>
                        </small>
                    );
                }

                if (isBeta) {
                    header.push(
                        <BetaTag />
                    );
                }

                return (
                    <Panel
                        header={header} collapsable={true}
                        expanded={this.state.collapsedItems[actualId]}
                        data-expander-id={actualId}
                        onSelect={onSelect}
                    >
                        { body }
                    </Panel>
                );
            })]);
        return (
            <Grid id="main-grid" className="help-page">
                {fragment === '' ? null :
                    <style>{`#${fragment} { animation-name: emphasis; animation-duration: 10s; } `}</style>}

                <Row ref={(me) => { if (me) { this.headerElem = $(me.getDOMNode()); } }} className="header-sticky">
                    <Col smOffset={1} sm={10} className="help-search-header">
                        {/*<h1>BRCA Exchange: Help</h1>*/}
                        <SearchController setExpansion={this.setExpansion} headerElem={this.headerElem} target="#help-body" />
                    </Col>
                </Row>

                <Row>
                    <Col smOffset={1} sm={10} id="help-body">
                        {helpTiles}
                    </Col>
                </Row>
            </Grid>
        );
    }
});

module.exports = Help;
