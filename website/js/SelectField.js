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
            <FormGroup controlId="formControlsSelect" className="d-flex align-items-center gap-2 mb-0">
                {label ? <FormLabel className="mb-0 text-nowrap fw-bold">{label}</FormLabel> : null}
                <FormSelect value={value} onChange={this.onChange}>
                    {opels}
                </FormSelect>
            </FormGroup>
        );
    }
}

export default SelectField;

