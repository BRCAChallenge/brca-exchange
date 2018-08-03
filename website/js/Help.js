/*eslint-env browser */
'use strict';

import classNames from "classnames";
const React = require('react');
const RawHTML = require('./RawHTML');
const {Grid, Col, Row, Panel, ListGroup, ListGroupItem, Glyphicon, CollapsableMixin, BootstrapMixin} = require('react-bootstrap');
const {State} = require('react-router');

const slugify = require('./slugify');
const content = require('./content');

const navbarHeight = 70;

const CollapsableListItem = React.createClass({
    mixins: [State, BootstrapMixin, CollapsableMixin],

    onClick(e) {
        this.setState({ expanded: !this.state.expanded });
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

    componentDidMount() {
        var fragment = slugify(window.location.hash.slice(1));
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

    render() {
        let help = localStorage.getItem("research-mode") === 'true' ? content.helpContentResearch : content.helpContentDefault;

        var fragment = slugify(window.location.hash.slice(1));
        // looks for the fragment within id attributes wrapped in quotes (in dev) or without quotes (minified)
        // we compile it once here and pass the same regex to each call to shouldBeExpanded() to save a little time
        const fragRegex = new RegExp(`(id=${fragment}|id="${fragment})"`);

        var helpTiles = help.map(({section, tiles}) =>
            [<h1>{section}</h1>, tiles.map(({name, id, contents, list, reference}) => {
                let header = [name];
                // if the user clicks a reference link in a tile header, don't toggle the tile, and open the link.
                let onSelect = function (e) {
                    if (e.target.classList.contains("help-reference-link")) {
                        e.selected = false;
                    } else {
                        e.preventDefault();
                    }
                };
                let body = [];
                if (contents) {
                    body.push(<RawHTML html={contents} />);
                }
                if (list) {
                    body.push(
                        <ListGroup fill>
                            { list.map(({name, id, contents}) =>
                                <CollapsableListItem defaultExpanded={this.shouldBeExpanded(fragment, fragRegex, {name, id, contents})}
                                                     id={id ? id : slugify(name)}
                                                     header={name}>
                                    <RawHTML html={contents} />
                                </CollapsableListItem>) }
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
                return (<Panel header={header} collapsable={true}
                               defaultExpanded={this.shouldBeExpanded(fragment, fragRegex, {name, id, contents, list})}
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
                    <Col smOffset={1} sm={10}>
                        {helpTiles}
                    </Col>
                </Row>
            </Grid>
        );
    }
});

module.exports = Help;
