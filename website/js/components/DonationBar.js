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
                    BRCA Exchange salutes National Hereditary Breast and Ovarian Cancer Week, 9/29/20 - 10/5/20.
                    To learn how you can help prevent heritable cancers, <Link to={`/donatedrive`}>click here</Link>.
	        	</p>
        	</div>
        );
    }
});

export default DonationBar;
