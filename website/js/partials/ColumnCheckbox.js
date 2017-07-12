/*global module: false, require: false */

import React from 'react';

import PureRenderMixin from '../helpers/PureRenderMixin'; // deep-equals version of PRM
import { Input } from 'react-bootstrap';

var ColumnCheckbox = React.createClass({
    mixins: [PureRenderMixin],
    onChange: function (e) {
        return this.props.onChange(e.target.value);
    },
    render: function () {
        var {label, title, initialCheck} = this.props;
        return (
            <div>
            <Input type="checkbox" label={title} checked={initialCheck[label]} onChange={this.onChange} />
            </div>
        );
    }
});

export default ColumnCheckbox;
