/*eslint-env browser */
/*global require: false */
'use strict';

import React from 'react';

const DonationBar = React.createClass({
    render() {
        return (
        	<div className="donation-bar">
    			<p>
                We are excited to introduce provisional assignments of the ACMG evidence codes for all variants in BRCA Exchange!  Please see our <a style={{textDecoration: 'underline'}} href="https://brcaexchange.org/blog/index.php/2024/08/14/new-feature-acmg-evidence-code-provisional-assignments/">latest blog post</a> for more information
	        	</p>
        	</div>
        );
    }
});

export default DonationBar;
