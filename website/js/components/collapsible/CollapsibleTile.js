/*eslint-env browser */
/*global require: false, module */
'use strict';

import React from "react";
import {Panel} from 'react-bootstrap';
import update from 'immutability-helper';

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
            if (React.isValidElement(c) && c.props.fieldName) {
                visibleFields[c.props.fieldName] = c.props.defaultVisible !== undefined ? c.props.defaultVisible : false;
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
        }), () => {
            // causes the parent to perform a (delayed) reflow
            this.props.onDimsChanged();
        });
    }

    fieldToggled(fieldName) {
        this.setState((pstate) => ({
            fieldExpansions: update(pstate.fieldExpansions, {
                [fieldName]: (visible) => visible ? !visible : true
            })
        }), () => {
            // causes the parent to perform a (delayed) reflow
            this.props.onDimsChanged();
        });
    }

    render() {
        const {fieldExpansions} = this.state;

        // TODO: figure out how to determine if everything is empty even though variant is in 10KG or ExAC
        const allEmpty = false;

        // create the source panel itself now
        const groupVisID = `group-panel-${this.props.groupTitle}`;
        const header = (
            <h3 style={{display: 'flex', flexDirection: 'row'}}>
                <a style={{flexGrow: 1}} href="#" onClick={(event) => this.props.onChangeGroupVisibility(groupVisID, event)}>
                    {this.props.groupTitle}
                </a>

                <a title='collapse all fields'
                    onClick={(event) => this.setAllFieldsExpansion(event, false)}
                    style={{cursor: 'pointer', marginRight: '10px'}}>
                    <i className="fa fa-angle-double-up" aria-hidden="true" />
                </a>

                <a title='expand all fields'
                    onClick={(event) => this.setAllFieldsExpansion(event, true)}
                    style={{cursor: 'pointer'}}>
                    <i className="fa fa-angle-double-down" aria-hidden="true" />
                </a>
            </h3>
        );

        // inject props for allowing us to react to child sections' request to change their visibility
        const togglableKids = React.Children.map(this.props.children, (c) => React.cloneElement(c, {
            onFieldToggled: this.fieldToggled,
            expanded: fieldExpansions[c.props.fieldName],
            hideEmptyItems: this.props.hideEmptyItems
        }));

        return (
            <div key={`group_collection-${groupVisID}`} className={ allEmpty && this.props.hideEmptyItems ? "group-empty variant-detail-group" : "variant-detail-group" }>
                <Panel
                    header={header}
                    collapsable={true}
                    defaultExpanded={localStorage.getItem("collapse-group_" + groupVisID) !== "true"}
                    hideEmptyItems={this.props.hideEmptyItems}>
                    { togglableKids }
                </Panel>
            </div>
        );
    }
};

CollapsibleTile.defaultProps = {
    onDimsChanged: () => {
        console.warn("onDimsChanged() unspecified; it should be specified to cause the parent container to reflow");
    }
};
