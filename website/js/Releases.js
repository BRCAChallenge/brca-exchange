'use strict';

var React = require('react');
var {Link} = require('react-router');
var {Table, Grid, Row, Col} = require('react-bootstrap');
var _ = require('underscore');
var backend = require('./backend');
var {dateFormat} = require('./util');
var config = require('config');

var Releases = React.createClass({
    getInitialState: () => ({ releases: {} }),
    componentWillMount: function() {
        backend.releases().subscribe(
            resp => this.setState(resp),
            () => this.setState({error: 'Problem connecting to server'}));
    },
    render: function () {
        /* eslint-disable dot-notation */
        var rows = _.map(this.state.releases, release => (
            <tr>
                <td><Link to={`/release/${release.id}`}>Version {release.id}</Link></td>
                <td>{dateFormat(release['date'])}</td>
                <td>{release['sources']}</td>
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
    getInitialState: () => ({ releases: [{}], latest: -1 }),
    componentWillMount: function() {
        backend.release(this.props.params.id).subscribe(
            resp => this.setState(resp),
            () => this.setState({error: 'Problem connecting to server'}));
    },
    componentWillReceiveProps: function(nextProps) {
        backend.release(nextProps.params.id).subscribe(
            resp => this.setState(resp),
            () => this.setState({error: 'Problem connecting to server'}));
    },
    render: function () {
        var release = this.state.releases[0],
            latest = this.state.latest,
            s = n => n === 1 ? '' : 's';
        if (release.hasOwnProperty('notes')) {
            var releaseNotes = release.notes.replace(/\n\s*\n/g, '\n\n');
        }
        /* eslint-disable dot-notation */
        return (
            <Grid fluid={true}>
                <Row>
                    <Col sm={8} smOffset={2} md={6} mdOffset={3} className='text-left'>
                        <h1>{release.id === latest && 'Current'} Release Notes</h1>
                        {release.id !== latest && <p>* Note that this is not the most current release. Click <Link to={`/release/${latest}`}>here</Link> to see the most current release.</p>}
                    </Col>
                </Row>
                <Row>
                    <Col sm={8} smOffset={2} md={6} mdOffset={3} className='text-center'>
                        <p className='release-notes text-left'>{releaseNotes}</p>
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
                                <td><a href={`${config.backend_url}/site_media/media/releases/${release['archive']}`}>Link</a></td>
                            </tr>
                            <tr>
                                <td className="active"><b>Date</b></td>
                                <td>{dateFormat(release['date'])}</td>
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
