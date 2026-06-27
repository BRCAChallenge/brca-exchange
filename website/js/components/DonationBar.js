/*eslint-env browser */
'use strict';

import React from 'react';

class DonationBar extends React.PureComponent {
    render() {
        return (
        	<div className="donation-bar">
    		    <p>
			Changes are coming to BRCA Exchange!  Please visit our <a style={{textDecoration: 'underline'}} href="https://brcaexchange.org/blog/index.php/2026/05/14/changes-are-coming-to-brca-exchange/">blog</a> for the latest news.
	        	</p>
        	</div>
        );
    }
}

export default DonationBar;
