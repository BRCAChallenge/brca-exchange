/*global module: false, require: false */
'use strict';

import {FormGroup, ControlLabel, FormControl} from 'react-bootstrap';

var React = require('react');
var PureRenderMixin = require('./PureRenderMixin'); // deep-equals version of PRM
var _ = require('underscore');

var SelectField = React.createClass({
    mixins: [PureRenderMixin],

    onChange: function (e) {
        return this.props.onChange(e.target.value);
    },

    render: function () {
        const {options, label, value} = this.props;
        const opels = _.map(options, v => <option key={v} value={v}>{v}</option>);

        return (
            <FormGroup controlId="formControlsSelect">
                <ControlLabel>{label}</ControlLabel>
                <FormControl value={value} componentClass="select" onChange={this.onChange} placeholder="select">
                    {opels}
                </FormControl>
            </FormGroup>
        );
    },
});

module.exports = SelectField;
