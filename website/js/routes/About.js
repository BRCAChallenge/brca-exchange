
var React = require('react');
var {Grid, Col, Row} = require('react-bootstrap');

var content = require('../data/content');
var RawHTML = require("../helpers/RawHTML");

var About = React.createClass({
    render: function() {
        var {page} = this.props.params;

        return (
            <Grid id="main-grid" className="main-grid">
                <Row>
                    <Col smOffset={1} sm={10}>
                        <RawHTML html={content.pages[page]} />
                    </Col>
                </Row>
            </Grid>
        );
    }
});

module.exports = About;
