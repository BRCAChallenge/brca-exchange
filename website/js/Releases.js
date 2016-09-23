/*global grecaptcha: false */
'use strict';

var React = require('react');
var {Grid, Row, Col} = require('react-bootstrap');
var _ = require('underscore');
var backend = require('./backend');

var Releases = React.createClass({
    getInitialState: () => ({ releases: {} }),
    componentWillMount: function() {
        backend.releases().subscribe(
            resp => this.setState({releases: resp}),
            () => this.setState({error: 'Problem connecting to server'}));
    },
    render: function () {
        var rows = _.map(this.state.releases, release => (
            <tr>
                <td>{release.id}</td>
                <td>{release.timestamp}</td>
                <td>{release.comment}</td>
                <td>{release.variants_added}</td>
                <td>meh</td>
                <td>{release.variants_modified}</td>
                <td>{release.variants_deleted}</td>
            </tr>
        ));
        return (
            <Grid fluid={true}>
                <Row>
                    <Col sm={8} smOffset={2}><h2>Release Notes</h2></Col>
                </Row>
                <Row>
                    <Col sm={8} smOffset={2}>
                        <table className="table table-bordered table-condensed">
                            <thead>
                                <th>Raw Data</th>
                                <th>Date</th>
                                <th>Data Sources</th>
                                <th>New Variants</th>
                                <th>New Classifications</th>
                                <th>Changed/Updated Variants</th>
                                <th>Deleted Variants</th>
                            </thead>
                            <tbody>
                                {rows}
                            </tbody>
                        </table>
                    </Col>
                </Row>
            </Grid>);
    }
});

var Release = React.createClass({
    getInitialState: () => ({}),
    render: function () {
        return (
            <Grid fluid={true}>
                <Row>
                    <Col sm={8} smOffset={2}>Release Notes</Col>
                </Row>
            </Grid>);
    }
});

module.exports = ({
    Releases: Releases,
    Release: Release,
});
