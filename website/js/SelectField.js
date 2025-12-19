/*global module: false, require: false */
'use strict';

import { FormGroup, FormLabel, FormSelect } from 'react-bootstrap';

import React from 'react';
var _ = require('underscore');

class SelectField extends React.PureComponent {

    onChange = (e) => this.props.onChange(e.target.value);

    render() {
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
    }
}

export default SelectField;

