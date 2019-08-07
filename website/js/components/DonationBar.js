/*eslint-env browser */
/*global require: false */
'use strict';

import React from 'react';

const DonationBar = React.createClass({
    render() {
        return (
        	<div className="donation-bar">
    			<p>
	        		Do you rely on high quality, open-source data from BRCA Exchange?
	        		Support our work by&nbsp;
	        		<a href="https://secure.ucsc.edu/s/1069/bp18/interior.aspx?sid=1069&gid=1001&pgid=780&cid=1749&dids=1004" target="_blank">
	        			donating today
	        		</a>!
	        	</p>
        	</div>
        );
    }
});

export default DonationBar;
