'use strict';

var React = require('react');
var {Link} = require('react-router');
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
                <td><a href={release.data_link}>Link</a></td>
                <td>{release.date_released}</td>
                <td>{release.data_sources}</td>
                <td><Link to={`/variants?release=${release.id}&changeTypes[]=new`}>{release.variants_added}</Link></td>
                <td><Link to={`/variants?release=${release.id}&changeTypes[]=added_classification&changeTypes[]=changed_classification`}>{release.variants_classified}</Link></td>
                <td><Link to={`/variants?release=${release.id}&changeTypes[]=added_information&changeTypes[]=changed_information`}>{release.variants_modified}</Link></td>
                <td><Link to={`/variants?release=${release.id}&changeTypes[]=deleted`}>{release.variants_deleted}</Link></td>
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
                                <th>Data</th>
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
