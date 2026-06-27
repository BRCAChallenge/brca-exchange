/*eslint-env browser */
'use strict';

import React from 'react';

class GroupHelpButton extends React.PureComponent {
    render() {
        const { onClick, style } = this.props;
        return (
            <span role='button' tabIndex={0} onClick={onClick} onKeyDown={(e) => {
		    if (e.key === 'Enter' || e.key === ' ') onClick?.(e);
	    }} aria-label="Help" style={style}
                className='panel-help-btn fa fa-question-circle'
            />
        );
    }
}

export default GroupHelpButton;
