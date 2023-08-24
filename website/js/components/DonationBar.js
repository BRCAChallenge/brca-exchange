/*eslint-env browser */
/*global require: false */
'use strict';

import React from 'react';

const DonationBar = React.createClass({
    render() {
        return (
        	<div className="donation-bar">
    			<p>
                Thank you to everyone who contributed to the BRCA Exchange user survey!  The results are available <a style={{textDecoration: 'underline'}} href="https://brcaexchange.org/blog/index.php/2023/08/20/results-of-the-brca-exchange-user-survey/">here</a> on the BRCA Exchange blog .
	        	</p>
        	</div>
        );
    }
});

export default DonationBar;
