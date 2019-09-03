/*eslint-env browser */
'use strict';

import classNames from "classnames";
import update from 'immutability-helper';
import BetaTag from "./components/BetaTag";
import {debounce} from "lodash";
// import debounce from 'lodash/debounce';
const React = require('react');
const RawHTML = require('./RawHTML');
const {Grid, Col, Row, Panel, ListGroup, ListGroupItem, Glyphicon, CollapsableMixin, BootstrapMixin} = require('react-bootstrap');
const {State} = require('react-router');

const slugify = require('./slugify');
const content = require('./content');

const navbarHeight = 70;

const RawHTMLHighlight = React.createClass({
    getInitialState() {
        return {
            matched: false,
            text: null
        };
    },

    simpleHighlight(query, orig) {
        if (!query) {
            return { text: orig, matches: 0 };
        }

        const r = new RegExp(query, "giu");
        let matches = 0;
        const text = orig.replace(r, match => {
            matches += 1;
            return `<span class="highlighted">${match}</span>`;
        });
        return { text, matches };
    },

    componentDidUpdate(prevProps, prevState) {
        if (this.props.searchTerm !== prevProps.searchTerm || this.props.html !== prevProps.html) {
            // recompute highlight, potentially expanding our parent if there's a match
            const result = this.simpleHighlight(this.props.searchTerm, this.props.html);
            this.setState({
                matched: result.matches > 0,
                text: result.text
            });
        }

        if (this.state.matched !== prevState.matched) {
            this.props.setExpansion(this.props.collapserId, this.state.matched);
        }
    },

    /*
    // also look up the collapser and trigger its click method
    this.props.setExpansion(this.props.collapserId, true);
     */

    render: function() {
        let content = this.state.matched ? this.state.text : this.props.html;

        return <RawHTML html={content} />;
    }
});
RawHTMLHighlight.propTypes = {
    html: React.PropTypes.string.isRequired,
    collapserId: React.PropTypes.string.isRequired,
    searchTerm: React.PropTypes.string
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
        this.debouncedCommitSearch = debounce(() => {
            console.log("Triggered");
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
                    console.log(e.target.getAttribute('id'), "clicked");

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
                    body.push(<RawHTMLHighlight html={contents} collapserId={actualId} setExpansion={this.setExpansion} searchTerm={this.state.committedSearch} />);
                }

                if (list) {
                    body.push(
                        <ListGroup fill>
                            {
                                list.map(({name, id, contents}) => {
                                    const actualId = id ? id : slugify(name);
                                    return (
                                        <CollapsableListItem
                                            id={actualId} setExpansion={this.setExpansion}
                                            expanded={this.state.collapsedItems[actualId]}
                                            header={name}>
                                            <RawHTMLHighlight html={contents} collapserId={actualId}  setExpansion={this.setExpansion} searchTerm={this.state.committedSearch} />
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

                return (<Panel header={header} collapsable={true}
                               expanded={this.state.collapsedItems[actualId]}
                               onSelect={onSelect}
                        >
                    { body }
                </Panel>);
            })]);
        return (
            <Grid id="main-grid" className="help-page">
                {fragment === '' ? null :
                    <style>{`#${fragment} { animation-name: emphasis; animation-duration: 10s; } `}</style>}

                <Row>
                  <Col smOffset={1} sm={10} className="help-search-header">
                    <h1>BRCA Exchange: Help</h1>
                    <div className="form-group has-feedback has-search">
                      <span className="glyphicon glyphicon-search form-control-feedback" />
                      <input type="text" className="form-control" placeholder="Search" value={this.state.searchTerm} onChange={this.searchChanged} />
                    </div>
                  </Col>
                </Row>

                <hr />

                <Row>
                    <Col smOffset={1} sm={10}>
                        {helpTiles}
                    </Col>
                </Row>
            </Grid>
        );
    }
});

module.exports = Help;
