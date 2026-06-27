/*eslint-env browser */
'use strict';

import React from "react";
import { Card, Collapse } from 'react-bootstrap';
import update from 'immutability-helper';

import GroupHelpButton from '../GroupHelpButton';

export default class CollapsibleTile extends React.Component {
    constructor(props) {
        super(props);

        // compute default group open state from localStorage (same logic as old defaultExpanded)
        const groupVisID = `group-panel-${props.groupTitle}`;
        const defaultExpanded = localStorage.getItem("collapse-group_" + groupVisID) !== "true";

        this.state = {
            fieldExpansions: this.getFieldExpansions(props),
            groupOpen: defaultExpanded
        };

        this.getFieldExpansions = this.getFieldExpansions.bind(this);
        this.fieldToggled = this.fieldToggled.bind(this);
        this.setAllFieldsExpansion = this.setAllFieldsExpansion.bind(this);
        this.toggleGroupOpen = this.toggleGroupOpen.bind(this);
    }

    getFieldExpansions(props) {
        if (!props.children) return {};

        const visibleFields = {};
        React.Children.map(props.children, (c) => {
            if (React.isValidElement(c) && c.props.id) {
                visibleFields[c.props.id] =
                    c.props.defaultVisible !== undefined ? c.props.defaultVisible : false;
            }
        });
        return visibleFields;
    }

    setAllFieldsExpansion(e, newExpansion) {
        e.stopPropagation();
        this.setState((pstate) => ({
            // generate a new fieldExpansions object where all existing fields' visibilities are set to newExpansion
            fieldExpansions: Object.keys(pstate.fieldExpansions)
                .reduce((o, k) => { o[k] = newExpansion; return o; }, {})
        }));
    }

    fieldToggled(id) {
        this.setState((pstate) => ({
            fieldExpansions: update(pstate.fieldExpansions, {
                [id]: (visible) => (visible ? !visible : true)
            })
        }));
    }

    toggleGroupOpen(event) {
        const groupVisID = `group-panel-${this.props.groupTitle}`;
        if (this.props.onChangeGroupVisibility) {
            this.props.onChangeGroupVisibility(groupVisID, event);
        }
        this.setState(prev => ({ groupOpen: !prev.groupOpen }));
    }

    render() {
        const { fieldExpansions } = this.state;

        // TODO: figure out how to determine if everything is empty even though variant is in 10KG or ExAC
        const allEmpty = this.props.allEmpty;

        // create the source panel itself now
        const groupVisID = `group-panel-${this.props.groupTitle}`;

        // inject props for allowing us to react to child sections' request to change their visibility
        const togglableKids = React.Children.map(this.props.children, (c) => {
            if (!React.isValidElement(c)) return c;

	    // Never inject custom props into DOM nodes like <div>
            if (typeof c.type === "string") return c;

            // Only inject into children that participate in field expansion (must have an id)
            const id = c.props && c.props.id;
            if (!id) return c;

	    return React.cloneElement(c, {
                onFieldToggled: this.fieldToggled,
                expanded: !!fieldExpansions[id],
                relayoutGrid: this.props.relayoutGrid,
                hideEmptyItems: this.props.hideEmptyItems
            });
	});

        return (
            <div
                key={`group_collection-${groupVisID}`}
                className={allEmpty && this.props.hideEmptyItems ? "group-empty variant-detail-group" : "variant-detail-group"}
            >
                <Card className="mb-3 shadow">
                    <Card.Header as="div" className="d-flex align-items-center fw-bold">
                        <span
                            role="button"
                            className="title text-decoration-none flex-grow-1"
                            onClick={this.toggleGroupOpen}
                            aria-expanded={this.state.groupOpen}
                        >
                            {this.props.displayTitle || this.props.groupTitle}
                        </span>

                        <span
                            title="collapse all fields"
                            className="toggle-subfields"
                            onClick={(event) => this.setAllFieldsExpansion(event, false)}
                            style={{ cursor: 'pointer', marginRight: '10px' }}
                        >
                            <i className="fa fa-angle-double-up" aria-hidden="true" />
                        </span>

                        <span
                            title="expand all fields"
                            className="toggle-subfields"
                            onClick={(event) => this.setAllFieldsExpansion(event, true)}
                            style={{ cursor: 'pointer' }}
                        >
                            <i className="fa fa-angle-double-down" aria-hidden="true" />
                        </span>

                        {this.props.helpSection && (
                            <GroupHelpButton
                                group={this.props.helpSection}
                                onClick={(event) => {
                                    this.props.showHelp(event, this.props.helpSection);
                                    return true;
                                }}
                            />
                        )}
                    </Card.Header>

                    <Collapse
                        in={this.state.groupOpen}
                        onEntered={this.props.relayoutGrid}
                        onExited={this.props.relayoutGrid}
                    >
                        <div>
                            <Card.Body>
                                {togglableKids}
                            </Card.Body>
                        </div>
                    </Collapse>
                </Card>
            </div>
        );
    }
};

CollapsibleTile.defaultProps = {
    relayoutGrid: () => {
        console.warn("relayoutGrid() unspecified; it should be specified to allow the parent container to reflow");
    }
};

