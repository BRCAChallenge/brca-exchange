/*eslint-env browser */
/*global require: false */
'use strict';

import React from 'react';

const DonationBar = React.createClass({
    render() {
        return (
        	<div className="donation-bar">
    			<p>
                We are committed to making BRCA Exchange more accessible to a general audience.  If you are a patient, previvor or other member of the general public, please visit our <a style={{textDecoration: 'underline'}}
	    href="https://docs.google.com/forms/d/e/1FAIpQLSec850bV5iIg3ZFxmPkcM4XeztLtoxQsjtgzqVNZXq_0KzjpA/viewform">survey!</a> Thank you for your time!
	        	</p>
        	</div>
        );
    }
});

export default DonationBar;
