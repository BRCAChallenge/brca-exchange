/*global module: false, require: false */

import React from 'react';

import PureRenderMixin from '../helpers/PureRenderMixin'; // deep-equals version of PRM
import { Input } from 'react-bootstrap';
import _ from 'underscore';

var SelectField = React.createClass({
    mixins: [PureRenderMixin],
    onChange: function (e) {
        return this.props.onChange(e.target.value);
    },
    render: function () {
        var {options, label, value} = this.props,
            opels = _.map(options, v => <option key={v} value={v}>{v}</option>);

        return (
            <Input type="select" value={value} label={label} onChange={this.onChange}>
                {opels}
            </Input>
        );
    }
});

export default SelectField;
