'use strict';

var React = require('react');
var PureRenderMixin = require('../helpers/PureRenderMixin'); // deep-equals version of PRM

var brcaLogo = require('../img/BRCA-Exchange-tall-tranparent.png');

var DisclaimerModal = require('./DisclaimerModal');

var Footer = React.createClass({
    mixins: [PureRenderMixin],
    render: function() {
        return (
            <div className="container footer">
                <div className="col-sm-5 left-footer">
                    <ul>
                        <li><a href="/">Home</a></li>
                        <li><a href="/about/history">About</a></li>
                        <li><a href="/variants">Variants</a></li>
                        <li><a href="/help">Help</a></li>
                        <li><a href="/about/api">API</a></li>
                    </ul>
                </div>
                <div className="col-sm-2 logo-footer">
                    <img href="#" src={brcaLogo} alt="brca exchange logo" />
                </div>
                <div className="col-sm-5 right-footer">
                    <ul>
                        <li>
                            <DisclaimerModal text="Disclaimer"/>
                        </li>
                        <li>
                            <a href="mailto:brca-exchange-contact@genomicsandhealth.org?subject=BRCA Exchange website">
                                Contact us
                            </a>
                        </li>
                        <li>
                            <a href="https://github.com/BD2KGenomics/brca-website">
                                Source code
                            </a>
                        </li>
                    </ul>
                </div>
            </div>
        );
    }
});

module.exports = Footer;
