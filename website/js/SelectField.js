/*global module: false, require: false */
'use strict';

import { FormGroup, FormLabel, FormSelect } from 'react-bootstrap';

var React = require('react');
var PureRenderMixin = require('./PureRenderMixin'); // deep-equals version of PRM
var _ = require('underscore');

var SelectField = React.createClass({
    mixins: [PureRenderMixin],

    onChange: function (e) {
        return this.props.onChange(e.target.value);
    },

    render: function () {
        const { options, label, value } = this.props;
        const opels = _.map(options, v => <option key={String(v)} value={v}>{v}</option>);

        return (
            <FormGroup controlId="formControlsSelect">
                {label ? <FormLabel>{label}</FormLabel> : null}
                <FormSelect value={value} onChange={this.onChange}>
                    {opels}
                </FormSelect>
            </FormGroup>
        );
    },
});

export default SelectField;

