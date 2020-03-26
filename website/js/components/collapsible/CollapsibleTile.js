/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";
import {Panel} from 'react-bootstrap';
import update from 'immutability-helper';

import GroupHelpButton from '../GroupHelpButton';

export default class CollapsibleTile extends React.Component {
    constructor(props) {
        super(props);

        this.state = {
            fieldExpansions: this.getFieldExpansions(props)
        };

        this.getFieldExpansions = this.getFieldExpansions.bind(this);
        this.fieldToggled = this.fieldToggled.bind(this);
        this.setAllFieldsExpansion = this.setAllFieldsExpansion.bind(this);
    }

    getFieldExpansions(props) {
        if (!props.children) {
            return {};
        }

        const visibleFields = {};

        React.Children.map(props.children, (c) => {
            if (React.isValidElement(c) && c.props.id) {
                visibleFields[c.props.id] = c.props.defaultVisible !== undefined ? c.props.defaultVisible : false;
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
                [id]: (visible) => visible ? !visible : true
            })
        }));
    }

    render() {
        const {fieldExpansions} = this.state;

        // TODO: figure out how to determine if everything is empty even though variant is in 10KG or ExAC
        const allEmpty = this.props.allEmpty;

        // create the source panel itself now
        const groupVisID = `group-panel-${this.props.groupTitle}`;

        // inject props for allowing us to react to child sections' request to change their visibility
        const togglableKids = React.Children.map(this.props.children, (c) => React.cloneElement(c, {
            onFieldToggled: this.fieldToggled,
            expanded: fieldExpansions[c.props.id],
            relayoutGrid: this.props.relayoutGrid,
            hideEmptyItems: this.props.hideEmptyItems
        }));

        return (
            <div key={`group_collection-${groupVisID}`} className={ allEmpty && this.props.hideEmptyItems ? "group-empty variant-detail-group" : "variant-detail-group" }>
                <Panel
                    collapsible={true}
                    defaultExpanded={localStorage.getItem("collapse-group_" + groupVisID) !== "true"}
                    hideEmptyItems={this.props.hideEmptyItems}
                >
                    <Panel.Heading>
                        <Panel.Title componentClass="h3">
                            <Panel.Toggle componentClass="a" className="title"
                                onClick={(event) => this.props.onChangeGroupVisibility(groupVisID, event)}
                            >
                                {this.props.displayTitle || this.props.groupTitle}
                            </Panel.Toggle>

                            <a title='collapse all fields'
                                className="toggle-subfields"
                                onClick={(event) => this.setAllFieldsExpansion(event, false)}
                                style={{cursor: 'pointer', marginRight: '10px'}}>
                                <i className="fa fa-angle-double-up" aria-hidden="true" />
                            </a>

                            <a title='expand all fields'
                                className="toggle-subfields"
                                onClick={(event) => this.setAllFieldsExpansion(event, true)}
                                style={{cursor: 'pointer'}}>
                                <i className="fa fa-angle-double-down" aria-hidden="true" />
                            </a>

                            {
                                this.props.helpSection &&
                                <GroupHelpButton group={this.props.helpSection}
                                    onClick={(event) => {
                                        this.props.showHelp(event, this.props.helpSection);
                                        return true;
                                    }}
                                />
                            }
                        </Panel.Title>
                    </Panel.Heading>
                    <Panel.Collapse
                        onEntered={this.props.relayoutGrid}
                        onExited={this.props.relayoutGrid}
                    >
                        <Panel.Body>
                        { togglableKids }
                        </Panel.Body>
                    </Panel.Collapse>

                </Panel>
            </div>
        );
    }
};

CollapsibleTile.defaultProps = {
    relayoutGrid: () => {
        console.warn("relayoutGrid() unspecified; it should be specified to allow the parent container to reflow");
    }
};
