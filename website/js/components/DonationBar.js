/*eslint-env browser */
/*global require: false */
'use strict';

import React from 'react';

const DonationBar = React.createClass({
    render() {
        return (
        	<div className="donation-bar">
    			<p>
                Please take our <a style={{textDecoration: 'underline'}} href="https://docs.google.com/forms/d/e/1FAIpQLScaFT5IxWP7zqCiCwAHMCIr-2CfvgyITQGaQS517c5S01lyfA/viewform">survey</a> on BRCA Exchange functionality (6 questions).
                    Your feedback is essential to help us improve the site and keep it freely available.
	        	</p>
        	</div>
        );
    }
});

export default DonationBar;
