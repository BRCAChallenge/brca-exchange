/*eslint-env browser */
'use strict';

import React from "react";
import update from 'immutability-helper';
import BetaTag from "./components/BetaTag";
import { debounce } from "lodash";
import { Container, Row, Col, Card, ListGroup, Collapse } from 'react-bootstrap';
import SearchController, { EXTRA_SEARCH_PADDING } from "./components/help/SearchController";
import HardlinkHelper, { idHelpClicked } from "./components/help/HardlinkHelper";

import RawHTML from './RawHTML'; // keep CJS require if RawHTML uses module.exports
const $ = require('jquery');
const slugify = require('./slugify');
const content = require('./content');

export const navbarHeight = 70;

/** Collapsible list item used within a help tile */
class CollapsableListItem extends React.Component {
  constructor(props) {
    super(props);
    this.onClick = this.onClick.bind(this);
  }

  onClick(e) {
    // handle inline id helper first
    if (idHelpClicked(e)) return;

    this.props.setExpansion(this.props.id);
    e.preventDefault();
  }

  render() {
    const { header, id, expanded, children, ...rest } = this.props;

    return (
      <ListGroup.Item {...rest}>
        <h4>
          <a
            href="#"
            onClick={this.onClick}
            style={{ color: "inherit", textDecoration: "none" }}
            id={id}
          >
            <small>
              <i
                className={`fa fa-chevron-${expanded ? 'down' : 'right'}`}
                aria-hidden="true"
              />{" "}
            </small>
            <span style={{ verticalAlign: "text-bottom" }}>{header}</span>
          </a>
          <HardlinkHelper id={id} />
        </h4>

        <Collapse in={expanded}>
          <div>
            {children}
          </div>
        </Collapse>
      </ListGroup.Item>
    );
  }
}

class Help extends React.Component {
  constructor(props) {
    super(props);

    const { fragment, fragRegex } = this.fragmentMatchers();
    const help = localStorage.getItem("research-mode") === 'true'
      ? content.helpContentResearch
      : content.helpContentDefault;

    // Build initial collapsed/expanded map
    const collapsedItems = {};
    help.forEach(({ tiles }) => {
      tiles.forEach(({ name, id, contents, list }) => {
        const key = id ? id : slugify(name);
        collapsedItems[key] = this.shouldBeExpanded(fragment, fragRegex, { name, id, contents, list });
        if (list) {
          list.forEach(({ name: lname, id: lid, contents: lcontents }) => {
            const lkey = lid ? lid : slugify(lname);
            collapsedItems[lkey] = this.shouldBeExpanded(fragment, fragRegex, { name: lname, id: lid, contents: lcontents });
          });
        }
      });
    });

    this.state = {
      searchTerm: '',
      collapsedItems
    };

    // Debounce search commits so we don't thrash the UI while typing
    this.debouncedCommitSearch = debounce(() => {
      this.setState(pstate => ({ committedSearch: pstate.searchTerm }));
    }, 300);

    this.setExpansion = this.setExpansion.bind(this);
    this.headerElem = null;
  }

  fragmentMatchers() {
    const fragment = slugify(window.location.hash.slice(1));
    return {
      fragment,
      // matches id=fragment or id="fragment" (minified/dev)
      fragRegex: new RegExp(`(id=${fragment}|id="${fragment}")`)
    };
  }

  componentDidMount() {
    const fragment = slugify(window.location.hash.slice(1));
    if (fragment !== '') {
      setTimeout(() => {
        const el = document.getElementById(fragment);
        if (el) {
          const elemTop = $(el).offset().top;
          const headerOffset = navbarHeight + $('.header-sticky').outerHeight() + EXTRA_SEARCH_PADDING;

          window.scrollTo({
            top: elemTop - headerOffset,
            behavior: 'smooth'
          });
        }
      }, 0);
    }
  }

  shouldBeExpanded(fragment, fragRegex, { name, id, contents, list }) {
    // If there is no fragment, don't auto-expand (except lists default to collapsed as before)
    if (!fragment) {
      return !!list;
    }

    const slug = id || slugify(name);
    if (slug === fragment) return true;
    if (contents && fragRegex.test(contents)) return true;
    if (list && list.some(elem => this.shouldBeExpanded(fragment, fragRegex, elem))) return true;
    return false;
  }

  setExpansion(id, forced = null) {
    this.setState(pstate => ({
      collapsedItems: update(pstate.collapsedItems, {
        [id]: (forced !== null) ? { $set: forced } : (visible) => (visible ? !visible : true)
      })
    }));
  }

  render() {
    const help = localStorage.getItem("research-mode") === 'true'
      ? content.helpContentResearch
      : content.helpContentDefault;

    const { fragment } = this.fragmentMatchers();

    const helpTiles = help.map(({ section, tiles }) => ([
      <h1 key={`section-${slugify(section)}`}>{section}</h1>,
      tiles.map(({ name, id, contents, list, reference, isBeta }) => {
        const actualId = id ? id : slugify(name);

        return (
          <Card key={actualId} className="mb-3" data-expander-id={actualId}>
            <Card.Header as="h4" className="d-flex align-items-center">
              <a
                className="identifier text-decoration-none"
                id={actualId}
                onClick={() => { this.setExpansion(actualId); }}
                role="button"
                aria-expanded={!!this.state.collapsedItems[actualId]}
              >
                {name}
              </a>

              <HardlinkHelper id={actualId} />

              {reference && (
                <small key="help_reference" className="ms-2">
                  <a href={reference} target="_blank" rel="noreferrer">
                    <i className="fa fa-link help-reference-link" aria-hidden="true" />
                  </a>
                </small>
              )}

              {isBeta && <BetaTag key="beta_tag" />}
            </Card.Header>

            <Collapse in={!!this.state.collapsedItems[actualId]}>
              <div>
                {contents && (
                  <Card.Body>
                    <RawHTML key="contents" hardlinks={true} html={contents} />
                  </Card.Body>
                )}

                {list && (
                  <ListGroup key="listgroup" variant="flush">
                    {list.map(({ name: lname, id: lid, contents: lcontents }) => {
                      const localId = lid ? lid : slugify(lname);
                      return (
                        <CollapsableListItem
                          key={localId}
                          id={localId}
                          setExpansion={this.setExpansion}
                          expanded={!!this.state.collapsedItems[localId]}
                          data-expander-id={localId}
                          header={lname}
                        >
                          <RawHTML hardlinks={true} html={lcontents} />
                        </CollapsableListItem>
                      );
                    })}
                  </ListGroup>
                )}
              </div>
            </Collapse>
          </Card>
        );
      })
    ]));

    return (
      <Container id="main-grid" className="help-page">
        {fragment === '' ? null : (
          <style>{`#${fragment} { animation-name: emphasis; animation-duration: 10s; }`}</style>
        )}

        <Row
          ref={(node) => {
            if (node) { this.headerElem = $(node); }
          }}
          className="header-sticky"
        >
          <Col sm={{ span: 10, offset: 1 }} className="help-search-header">
            <SearchController
              researchMode={localStorage.getItem('research-mode')}
              setExpansion={this.setExpansion}
              headerElem={this.headerElem}
              target="#help-body"
            />
          </Col>
        </Row>

        <Row>
          <Col sm={{ span: 10, offset: 1 }} id="help-body">
            {helpTiles}
          </Col>
        </Row>
      </Container>
    );
  }
}

export default Help;

