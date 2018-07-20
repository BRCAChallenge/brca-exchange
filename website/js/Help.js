/*eslint-env browser */
'use strict';

const React = require('react');
const RawHTML = require('./RawHTML');
const {Grid, Col, Row, Panel, ListGroup, ListGroupItem, Glyphicon} = require('react-bootstrap');
const {State} = require('react-router');

const slugify = require('./slugify');
const content = require('./content');

const navbarHeight = 70;

const CollapsibleListItem = React.createClass({
    mixins: [State],

    getInitialState() {
        return {
            expanded: !!this.props.defaultExpanded
        };
    },

    onClick(e) {
        this.setState({ expanded: !this.state.expanded });
        e.preventDefault();
    },

    render: function() {
        let {header, id, ...rest} = this.props;
        header = (
            <h4>
                <a href="#" onClick={this.onClick}
                   style={{color: "inherit", textDecoration: "none"}}
                   id={id}>{header}</a>
            </h4>);
        return (
            <ListGroupItem header={header} {...rest}>
                <div className={this.state.expanded ? "" : "hidden"} >
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

    shouldBeExpanded(fragment, {name, id, contents, list}) {
        let slug = id || slugify(name);
        if (slug === fragment) {
            return true;
        }
        if (contents && contents.includes(`id="${fragment}"`)) {
            return true;
        }
        if (list && list.some(elem => this.shouldBeExpanded(fragment, elem))) {
            return true;
        }
        return false;
    },

    render() {
        var fragment = slugify(window.location.hash.slice(1));

        var helpTiles = content.helpContent.map(({section, tiles}) =>
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
                                <CollapsibleListItem defaultExpanded={this.shouldBeExpanded(fragment, {name, id, contents})}
                                                     id={id ? id : slugify(name)}
                                                     header={name}>
                                    <RawHTML html={contents} />
                                </CollapsibleListItem>) }
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
                               defaultExpanded={this.shouldBeExpanded(fragment, {name, id, contents, list})}
                               onSelect={onSelect}
                        >
                    { body }
                </Panel>);
            })]);
        return (
            <Grid id="main-grid" className="help">
                {fragment === '' ? null :
                    <style>{`#${fragment} { animation-name: emphasis; animation-duration: 10s; } `}</style>}
                <Row>
                    <Col smOffset={1} sm={10}>
                        {helpTiles}
                    </Col>
                </Row>
            </Grid>
        );
    /*
        var helpContent;
        if (localStorage.getItem("research-mode") === 'true') {
            helpContent = content.pages.helpResearch;
        } else {
            helpContent = content.pages.help;
        }
        return (
            <Grid id="main-grid" className="help">
                {fragment === '' ? null :
                    <style>{`#${fragment} { animation-name: emphasis; animation-duration: 10s; } `}</style>}
                <Row>
                    <Col smOffset={1} sm={10}>
                        <RawHTML ref='content' html={helpContent}/>
                    </Col>
                </Row>
            </Grid>
        );
    */
    }
});

module.exports = Help;
