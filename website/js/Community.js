'use strict';

var React = require('react');
var PureRenderMixin = require('./PureRenderMixin'); // deep-equals version of PRM
var {Grid, Col, Row, Button, Table} = require('react-bootstrap');
var backend = require('backend');
var config  = require('./config');
var {Navigation, Link} = require('react-router');
var {Pagination} = require('react-data-components-bd2k');
var _ = require('underscore');
var placeholder = require('./img/placeholder.png');
var auth = require('./auth');

var Community = React.createClass({
    mixins: [PureRenderMixin, Navigation],
    componentWillMount: function () {
        this.fetch(this.state);
    },
    getInitialState: function () {
        return {
            data: [],
            pageLength: 10,
            page: 0,
            totalPages: 1,
            windowWidth: window.innerWidth
        };
    },
    setPages: function({data, count}) {
        return {
            data,
            count,
            totalPages: Math.ceil(count / this.state.pageLength)
        };
    },
    fetch: function (state) {
        backend.users(state).subscribe(
            resp => this.setState(this.setPages(resp)), // set data, count, totalPages
            () => this.setState({error: 'Problem connecting to server'}));
    },
    setStateFetch: function (opts) {
        var newState = {...this.state, ...opts};
        this.setState(newState);
        this.fetch(newState);
    },
    onChangePage: function (pageNumber) {
        this.setStateFetch({page: pageNumber});
    },
    logout: function() {
        auth.logout();
        this.forceUpdate();
    },
    render: function () {
        var queryParams = this.context.router.getCurrentQuery();
        var message;
        if (queryParams.registrationSuccess) {
            message =  (
				<div className="alert alert-success">
					<p>Thanks for signing up. We have sent you an email with a confirmation link to complete your registration.</p>
				</div>);
        }
        var {data, page, totalPages, error} = this.state;
        var rows = _.map(data, row => {

            var avatar;
            if (row.has_image) {
                var avatarLink = config.backend_url + '/site_media/media/' + row.id;
                avatar = <object className="avatar" data={avatarLink} type="image/jpg"/>;
            } else {
                avatar = <img src={placeholder}/>;
            }

            var {city, state, country} = row;
            var locationString = _.values(_.pick({city, state, country}, v => v)).join(', ');

			/*eslint-disable dot-notation*/
            return (
				<tr>
					<td>
						{avatar}
					</td>
					<td>
						<span id="name" className="row-wrap"><h3>{row['firstName']} {row['lastName']}, {row['title']}</h3></span>
						<span id="affiliation" className="row-wrap"><h4>{row['affiliation']}</h4></span>
						<span id="institution" className="row-wrap">{row['institution']}</span>
						<span id="location" className="row-wrap">{locationString}</span>
						<span id="contact" className="row-wrap">{row['email']} {row['phone_number']}</span>
					</td>
				</tr>);
			/*eslint-enable dot-notation*/
        });

        return (error ? <p>{error}</p> :
            <Grid id="main-grid">
                <Row id="message"> {message} </Row>
                <Row>
                    <Col smOffset="5">
                        <h3>BRCA Community</h3>
                    </Col>
                    {!auth.loggedIn() && <Col sm={1} smOffset="2"><Link to="/signup"><Button bsStyle="link">Sign up </Button></Link></Col>}
                    {!auth.loggedIn() && <Col sm={1}><Link to="/signin"><Button bsStyle="link"> Sign in </Button></Link></Col>}

                    {auth.loggedIn() && <Col sm={1} smOffset="2"><Link to="/profile"><Button bsStyle="link">Edit profile</Button></Link></Col>}
                    {auth.loggedIn() && <Col sm={1}><Button onClick={this.logout} bsStyle="link">Sign out</Button></Col>}

                </Row>
                <Row className="btm-buffer">
                    <Col sm={10}>
                        <Pagination
                            className="pagination pull-right-sm"
                            currentPage={page}
                            totalPages={totalPages}
                            onChangePage={this.onChangePage} />
                    </Col>
                </Row>
                <Row>
                    <Col md={8} mdOffset={2}>
                        <Table striped bordered>
                            <tbody>
                                {rows}
                            </tbody>
                        </Table>
                    </Col>
                </Row>
                <Row/>
            </Grid>
        );
    }
});

module.exports = Community;
