/*eslint-env browser */
'use strict';

const React = require('react');
const {State} = require('react-router');
const {Grid, Row, Col} = require('react-bootstrap');
const content = require('./content');

const WhyDonate = React.createClass({
    mixins: [State],

    render() {
        return (
            <Grid id="main-grid" className="help-page">
                <Row>
                    <Col smOffset={1} sm={10}>
                        {content.whyDonate}
                    </Col>
                </Row>
            </Grid>
        );
    }
});

module.exports = WhyDonate;
