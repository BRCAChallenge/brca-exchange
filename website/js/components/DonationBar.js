/*eslint-env browser */
/*global require: false */
'use strict';

import React from 'react';

const DonationBar = React.createClass({
    render() {
        return (
        	<div className="donation-bar">
    			<p>
                We are pleased to introduce a new feature: persistent URLs for the variants in BRCA Exchange!  More information is available <a style={{textDecoration: 'underline'}} href="https://brcaexchange.org/blog/index.php/2024/02/01/new-feature-persistent-urls-for-brca-exchange-variants/">here</a> on the BRCA Exchange blog .
	        	</p>
        	</div>
        );
    }
});

export default DonationBar;
