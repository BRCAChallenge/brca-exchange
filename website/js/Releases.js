'use strict';

var React = require('react');
var {Link} = require('react-router');
var {Table, Grid, Row, Col} = require('react-bootstrap');
var _ = require('underscore');
var backend = require('./backend');
var {dateFormat} = require('./util');

var Releases = React.createClass({
    getInitialState: () => ({ releases: {} }),
    componentWillMount: function() {
        backend.releases().subscribe(
            resp => this.setState({releases: resp}),
            () => this.setState({error: 'Problem connecting to server'}));
    },
    render: function () {
        /* eslint-disable dot-notation */
        var rows = _.map(this.state.releases, release => (
            <tr>
                <td><Link to={`/release/${release.id}`}>Version {release.id}</Link></td>
                <td>{dateFormat(release['date_released'])}</td>
                <td>{release['data_sources']}</td>
                <td><Link to={`/variants?release=${release.id}&changeTypes[]=new`}>{release['variants_added']}</Link></td>
                <td><Link to={`/variants?release=${release.id}&changeTypes[]=added_classification&changeTypes[]=changed_classification`}>{release['variants_classified']}</Link></td>
                <td><Link to={`/variants?release=${release.id}&changeTypes[]=added_information&changeTypes[]=changed_information`}>{release['variants_modified']}</Link></td>
                <td><Link to={`/variants?release=${release.id}&changeTypes[]=deleted`}>{release['variants_deleted']}</Link></td>
            </tr>
        ));
        /* eslint-enable dot-notation */
        return (
            <Grid fluid={true}>
                <Row>
                    <Col sm={8} smOffset={2}><h2>Release Notes</h2></Col>
                </Row>
                <Row>
                    <Col sm={8} smOffset={2}>
                        <table className="table table-bordered table-condensed">
                            <thead>
                                <th>Notes</th>
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
    getInitialState: () => ({release: {}}),
    componentWillMount: function() {
        backend.release(this.props.params.id).subscribe(
            resp => this.setState({release: resp[0]}),
            () => this.setState({error: 'Problem connecting to server'}));
    },
    render: function () {
        var release = this.state.release;
        var s = n => n === 1 ? '' : 's';
        /* eslint-disable dot-notation */
        return (
            <Grid fluid={true}>
                <Row>
                    <Col sm={8} smOffset={2} className='text-center'>
                        <h1>{release['is_current'] && 'Current'} Release Notes</h1>
                        <span>{release['release_notes']}</span>
                        <h3><Link to={`/variants?release=${release.id}&changeTypes[]=new`}>{release['variants_added']} new variant{s(release['variants_added'])}</Link></h3>
                        <h3><Link to={`/variants?release=${release.id}&changeTypes[]=added_classification&changeTypes[]=changed_classification`}>{release['variants_classified']} new classification{s(release['variants_classified'])}</Link></h3>
                        <h3><Link to={`/variants?release=${release.id}&changeTypes[]=added_information&changeTypes[]=changed_information`}>{release['variants_modified']} changed/updated variant{s(release['variants_modified'])}</Link></h3>
                        <h3><Link to={`/variants?release=${release.id}&changeTypes[]=deleted`}>{release['variants_deleted']} deleted variant{s(release['variants_deleted'])}</Link></h3>
                    </Col>
                </Row>
                <Row>
                    <Col sm={4} smOffset={4}>
                        <Table bordered>
                            <tr>
                                <td className="active"><b>Link to Data</b></td>
                                <td><a href={release['data_link']}>Link</a></td>
                            </tr>
                            <tr>
                                <td className="active"><b>md5sum of data</b></td>
                                <td>{release['md5sum']}</td>
                            </tr>
                            <tr>
                                <td className="active"><b>Date</b></td>
                                <td>{dateFormat(release['date_released'])}</td>
                            </tr>
                            <tr>
                                <td className="active"><b>Link to Data Schema</b></td>
                                <td><a href='#'>Link</a></td>
                            </tr>
                        </Table>
                {release['data_sources']}
                    </Col>
                </Row>
            </Grid>);
            /* eslint-enable dot-notation */
    }
});

module.exports = ({
    Releases: Releases,
    Release: Release,
});
