/*eslint-env browser */
'use strict';

import classNames from "classnames";
import update from 'immutability-helper';
import BetaTag from "./components/BetaTag";
import {debounce} from "lodash";
import SearchController, {EXTRA_SEARCH_PADDING} from "./components/help/SearchController";
import HardlinkHelper, {idHelpClicked} from "./components/help/HardlinkHelper";
// import debounce from 'lodash/debounce';
const React = require('react');
const ReactDOM = require('react-dom');
const RawHTML = require('./RawHTML');
const {Grid, Col, Row, Panel, ListGroup, ListGroupItem, Glyphicon, CollapsableMixin, BootstrapMixin} = require('react-bootstrap');
const {State} = require('react-router');

var $ = require('jquery');
const slugify = require('./slugify');
const content = require('./content');

export const navbarHeight = 70;

const CollapsableListItem = React.createClass({
    mixins: [State, BootstrapMixin, CollapsableMixin],

    onClick(e) {
        // handle the inline id helper if that's what the user clicked
        if (idHelpClicked(e)) { return; }

        // toggle our expansion status via this callback declared in our parent
        this.props.setExpansion(this.props.id);
        e.preventDefault();
    },

    getCollapsableDimensionValue() {
        return ReactDOM.findDOMNode(this.refs.content).scrollHeight;
    },

    getCollapsableDOMNode() {
        if (!this.isMounted() || !this.refs || !this.refs.content) {
            return null;
        }

        return ReactDOM.findDOMNode(this.refs.content);
    },

    render: function() {
        let {header, id, ...rest} = this.props;

        const headerElem = (
            <h4>
                <a href="#" onClick={this.onClick}
                    style={{color: "inherit", textDecoration: "none"}}
                    id={id} >
                    <small><Glyphicon  glyph={this.props.expanded ? "chevron-down" : "chevron-right"} /> </small>
                    <span style={{verticalAlign: "text-bottom"}}>{header}</span>
                </a>
                <HardlinkHelper id={id} />
            </h4>
        );

        return (
            <ListGroupItem header={headerElem} {...rest}>
                <div className={classNames(this.getCollapsableClassSet("collapse"))} ref="content">
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
                    const elemTop = $(el).offset().top;
                    const headerOffset = navbarHeight + $('.header-sticky').outerHeight() + EXTRA_SEARCH_PADDING;
                    console.log("Element top: ", elemTop, "; header offset: ", headerOffset, "; total: ", elemTop - headerOffset);
                    window.scrollTo({
                        top: elemTop - headerOffset,
                        behavior: 'smooth'
                    });
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

    onSelect(e) {
        // handle the inline id helper if that's what the user clicked
        if (idHelpClicked(e)) { return; }

        // find the parent that has an identifier associated with it
        const targetID = e.target.classList.contains("identifier")
            ? e.target.getAttribute('id')
            : $(e.target).parent('.identifier').attr('id');
        this.setExpansion(targetID);

        if (e.target.classList.contains("help-reference-link")) {
            // if the user clicks a reference link in a tile header, don't toggle the tile, and open the link.
            e.selected = false;
        } else {
            e.preventDefault();
        }
    },

    render() {
        let help = localStorage.getItem("research-mode") === 'true' ? content.helpContentResearch : content.helpContentDefault;

        const {fragment} = this.fragmentMatchers();

        var helpTiles = help.map(({section, tiles}) =>
            [<h1>{section}</h1>, tiles.map(({name, id, contents, list, reference, isBeta}) => {
                const actualId = id ? id : slugify(name);

                let header = [<span key="header_name" className="identifier" id={actualId}>{name}</span>];
                header.push(<HardlinkHelper id={actualId} />);

                let body = [];
                if (contents) {
                    body.push(<RawHTML key="contents" hardlinks={true} html={contents} />);
                }

                if (list) {
                    body.push(
                        <ListGroup key="listgroup" fill>
                            {
                                list.map(({name, id, contents}) => {
                                    const localId = id ? id : slugify(name);
                                    return (
                                        <CollapsableListItem
                                            key={localId} id={localId} setExpansion={this.setExpansion}
                                            expanded={this.state.collapsedItems[localId]}
                                            data-expander-id={localId}
                                            header={name}>
                                            <RawHTML hardlinks={true} html={contents} />
                                        </CollapsableListItem>
                                    );
                                })
                            }
                        </ListGroup>
                    );
                }

                if (reference) {
                    header.push(
                        <small key="help_reference">&nbsp;
                            <a href={reference} target="_blank">
                                <Glyphicon glyph="link" className="help-reference-link" />
                            </a>
                        </small>
                    );
                }

                if (isBeta) {
                    header.push(
                        <BetaTag key="beta_tag" />
                    );
                }

                return (
                    <Panel
                        header={header} collapsable={true}
                        key={actualId}
                        expanded={this.state.collapsedItems[actualId]}
                        data-expander-id={actualId}
                        onSelect={this.onSelect}
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
                        <SearchController
                            researchMode={localStorage.getItem('research-mode')}
                            setExpansion={this.setExpansion}
                            headerElem={this.headerElem} target="#help-body"
                        />
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
