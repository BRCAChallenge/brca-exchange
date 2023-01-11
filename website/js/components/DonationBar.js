/*eslint-env browser */
/*global require: false */
'use strict';

import React from 'react';

const {Link} = require('react-router');

const DonationBar = React.createClass({
    render() {
        return (
        	<div className="donation-bar">
    			<p>
                    Please take our <a href="https://docs.google.com/forms/d/e/1FAIpQLScaFT5IxWP7zqCiCwAHMCIr-2CfvgyITQGaQS517c5S01lyfA/viewform">brief survey</a> on BRCA Exchange usage.
                    Your feedback will help us improve the site and keep it freely available.
	        	</p>
        	</div>
        );
    }
});

export default DonationBar;
