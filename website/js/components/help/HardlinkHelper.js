'use strict';

import React from 'react';

export const HARDLINK_GLYPH = `Copy URL`;

/**
 * Handles clicking link helpers, which should change the address bar without performing an actual navigation.
 *
 * Note that this method is used not only in this component, but by the ad-hoc hardlink helpers created in RawHTML.js.
 * @param e the click event object
 * @returns {boolean} true if handled here, false otherwise
 */
export function idHelpClicked(e) {
    if (e.target.classList.contains("id-helper-tip")) {
        history.pushState(null, null, e.target.getAttribute('href'));

        // create a temporary DOM element into which we'll do a copy of the current URL
        const dummy = document.createElement('input'), text = window.location.href;
        document.body.appendChild(dummy);
        dummy.value = text;
        dummy.select();
        document.execCommand('copy');
        document.body.removeChild(dummy);

        e.preventDefault();
        return true;
    }

    return false;
}

const HardlinkHelper = React.createClass({
    render() {
        const {id} = this.props;

        return (
            <a className="id-helper-tip" onClick={(e) => { idHelpClicked(e); }} href={`#${id}`}
                dangerouslySetInnerHTML={{__html: HARDLINK_GLYPH}}
            />
        );
    }
});

export default HardlinkHelper;
