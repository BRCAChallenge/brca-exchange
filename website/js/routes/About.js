
import React from 'react';
import { Grid, Col, Row } from 'react-bootstrap';
import content from '../data/content';
import RawHTML from '../helpers/RawHTML';

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

export default About;
