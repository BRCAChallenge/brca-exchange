import React from 'react';
import { Navigation } from 'react-router';
import { Grid, Col, Row, Button, Modal, Glyphicon } from 'react-bootstrap';
import _ from 'underscore';
import content from '../data/content';
import logos from '../data/logos';
import VariantSearch from '../partials/VariantSearch';
import RawHTML from '../helpers/RawHTML';

var Home = React.createClass({
    mixins: [Navigation],
    getInitialState() {
        return {
            index: 0,
            direction: null,
            showModal: false
        };
    },

    onSearch(value) {
        this.transitionTo('/variants', null, {search: value});
    },
    render: function() {
        var {suggestions} = this.props;
        var logoItems = _.map(logos, ({id, logo, url}) => (
            <Col key={id} lg={4} md={6} xs={12} className="logo-item">
                <a href={url}>
                    <img id={id} src={logo} alt={id + ' logo'} />
                </a>
            </Col>
        ));
        return (
            <Grid id="main-grid" className='home'>
                <Row>
                    <Col smOffset={2} sm={8}>
                        <VariantSearch
                            id='home-search'
                            suggestions={suggestions}
                            onSearch={this.onSearch}/>
                    </Col>
                </Row>
                <Row>
                    <div className="jumbotron">
                        <RawHTML html={content.pages.home} />
                        <Button bsStyle="primary" className="center-block video-button" onClick={() => this.setState({ showModal: true })}>
                            <Glyphicon glyph="play-circle" />&nbsp;&nbsp;Video Overview
                        </Button>
                    </div>
                </Row>
                <Row className="logo-block">
                    {logoItems}
                </Row>
                {this.state.showModal && <Modal bsSize="large" onRequestHide={() => this.setState({ showModal: false })}>
                    <iframe className="vimeo-video" src="https://player.vimeo.com/video/199396428" frameBorder="0" webkitallowfullscreen mozallowfullscreen allowFullScreen>
                    FIXME: no iframe support
                    </iframe>
                </Modal>}
            </Grid>
        );
    }
});

export default Home;
