'use strict';

var _ = require('underscore');
var moment = require('moment');

var React = require('react');
var {Link} = require('react-router');
var {Table, Grid, Row, Col} = require('react-bootstrap');

var backend = require('./backend');
var config = require('config');
var anchorme = require("anchorme").default;


var Releases = React.createClass({
    getInitialState: () => ({ releases: {} }),
    componentWillMount: function() {
        backend.releases().subscribe(
            resp => this.setState(resp),
            () => this.setState({error: 'Problem connecting to server'}));
    },
    getSourceRepresentations: (sources) => {
        // exLOVD was renamed ExUV in October 2017
        return sources.replace(/exlovd/ig, 'ExUV');
    },
    render: function () {
        // Ensure releases are in descending order
        var releases = this.state.releases;
        if (Array.isArray(releases)) {
            releases = releases.sort(function(a, b) {
                return a.name - b.name;
            }).reverse();
        }
        var rows = _.map(releases, release => (
            <tr key={release.id}>
                <td style={{ whiteSpace: 'nowrap' }}><Link to={`/release/${release.id}`}>Version {release.name}</Link></td>
                <td style={{ whiteSpace: 'nowrap' }}>{moment(release.date, "YYYY-MM-DDTHH:mm:ss").format("DD MMMM YYYY")}</td>
                <td>{this.getSourceRepresentations(release.sources)}</td>
                <td><Link to={`/variants?release=${release.id}&changeTypes[]=new`}>{release['variants_added']}</Link></td>
                <td><Link to={`/variants?release=${release.id}&changeTypes[]=added_classification&changeTypes[]=changed_classification`}>{release['variants_classified']}</Link></td>
                <td><Link to={`/variants?release=${release.id}&changeTypes[]=added_information&changeTypes[]=changed_information`}>{release['variants_modified']}</Link></td>
                <td><Link to={`/variants?release=${release.id}&changeTypes[]=deleted`}>{release['variants_deleted']}</Link></td>
            </tr>
        ));
        return (
            <Grid id="main-grid" className="main-grid">
                <Row>
                    <Col smOffset={1} sm={10}>
                        <h1>Release Notes</h1>
                        <br />
                        <div className="table-responsive table-responsive-outset">
                            <table className="table table-inset-bordered table-grayheader nopointer">
                                <thead>
                                    <tr>
                                        <th>Notes</th>
                                        <th>Date</th>
                                        <th>Data Sources</th>
                                        <th>New Variants</th>
                                        <th>New Classifications</th>
                                        <th>Changed/Updated Variants</th>
                                        <th>Removed Variants</th>
                                    </tr>
                                </thead>
                                <tbody>
                                    {rows}
                                </tbody>
                            </table>
                        </div>
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
    generateReleaseNotes: function() {
        var release = this.state.releases[0];
        var releaseNotes = '';
        if (release.hasOwnProperty('notes')) {
            // format linebreaks
            releaseNotes = release.notes.replace(/\n\s*\n/g, '\n\n');
            // format hyperlinks
            releaseNotes = anchorme(releaseNotes);
            // rename exLOVD ExUV following the name change in October 2017
            releaseNotes = releaseNotes.replace(/exlovd/ig, 'ExUV');
        }
        return {__html: releaseNotes};
    },
    render: function () {
        var release = this.state.releases[0],
            latest = this.state.latest,
            s = n => n === 1 ? '' : 's';
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
                        <p className='release-notes text-left' dangerouslySetInnerHTML={this.generateReleaseNotes()}></p>
                        <h3><Link to={`/variants?release=${release.id}&changeTypes[]=new`}>{release['variants_added']} new variant{s(release['variants_added'])}</Link></h3>
                        <h3><Link to={`/variants?release=${release.id}&changeTypes[]=added_classification&changeTypes[]=changed_classification`}>{release['variants_classified']} new classification{s(release['variants_classified'])}</Link></h3>
                        <h3><Link to={`/variants?release=${release.id}&changeTypes[]=added_information&changeTypes[]=changed_information`}>{release['variants_modified']} changed/updated variant{s(release['variants_modified'])}</Link></h3>
                        <h3><Link to={`/variants?release=${release.id}&changeTypes[]=deleted`}>{release['variants_deleted']} removed variant{s(release['variants_deleted'])}</Link></h3>
                    </Col>
                </Row>
                <Row>
                    <Col sm={4} smOffset={4}>
                        <Table bordered>
                            <tr>
                                <td className="active"><b>Link to Data</b></td>
                                <td><a href={release.archive ? `${config.backend_url}/downloads/releases/${release.archive.split('.')[0]}/${release.archive}` : ''}>Link</a></td>
                            </tr>
                            <tr>
                                <td className="active"><b>Date</b></td>
                                <td>{moment(release.date, "YYYY-MM-DDTHH:mm:ss").format("DD MMMM YYYY")}</td>
                            </tr>
                        </Table>
                {release['data_sources']}
                    </Col>
                </Row>
            </Grid>);
    }
});

module.exports = ({
    Releases: Releases,
    Release: Release,
});
