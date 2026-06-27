'use strict';

import React from 'react';
import { Form } from 'react-bootstrap';

class ColumnCheckbox extends React.PureComponent {
    render() {
        const {label, title, initialCheck, onChange} = this.props;
        return (
            <div>
                <Form.Check type="checkbox" label={title} checked={!!(initialCheck && initialCheck[label])} onChange={onChange}/>
            </div>
        );
    }
}

export default ColumnCheckbox;
