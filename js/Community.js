'use strict';

var React = require('react');
var PureRenderMixin = require('./PureRenderMixin'); // deep-equals version of PRM
var {Grid, Col, Row, Button, Table} = require('react-bootstrap');
var backend = require('backend');
var config  = require('./config')
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
            message =  <div className="alert alert-success">
                <p>Thanks for signing up. We have sent you an email with a confirmation link to complete your registration.</p>
            </div>
        }
        var {data, page, totalPages, error} = this.state;
        var rows = _.map(data, row => {

            var avatar;
            if (row.has_image) {
                var avatar_link = config.backend_url + '/site_media/media/' + row['id']
                avatar = <object className="avatar" data={avatar_link} type="image/jpg"/>
            } else {
                avatar = <img src={placeholder}/>
            }

            var {city, state, country} = row;
            var location_string = _.values(_.pick({city,state,country}, v => v)).join(', ')

            return <tr >
                <td>
                    {avatar}
                </td>
                <td>
                    <span id="name" className="row-wrap"><h3>{row['firstName']} {row['lastName']}, {row['title']}</h3></span>
                    <span id="affiliation" className="row-wrap"><h4>{row['affiliation']}</h4></span>
                    <span id="institution" className="row-wrap">{row['institution']}</span>
                    <span id="location" className="row-wrap">{location_string}</span>
                    <span id="contact" className="row-wrap">{row['email']} {row['phone_number']}</span>
                </td>
            </tr>
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
                <Row>
                    <Col md={6} mdOffset={3} style={{height: "486px"}}>
                        <CommunityMap />
                    </Col>
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
                <Row>
                </Row>
            </Grid>
        );
    }
});

var CommunityMap = React.createClass({
    loadMap: function() {
        var map;
        window.initMap = function () {
            this.geo = new google.maps.Geocoder();
            map = window.map = new google.maps.Map(document.getElementById('communityMap'), {
                center: {lat: -0, lng: 0},
                zoom: 1
            });
            backend.userLocations().subscribe(({data}) => {
                _.map(data, ({firstName, lastName, title, institution, city, state, country})  => {
                    this.geo.geocode({address: city + "," + state + "," + country}, (results, status) => {
                        var l = results[0].geometry.location;
                        var userInfo = <div>{firstName} {lastName} {title}<br/>{institution}</div>;
                        var marker = new google.maps.Marker({position: { lat: l.lat(), lng: l.lng() }, map: map, title: "TEST LOCATION"});
                        var info = new google.maps.InfoWindow({content: React.renderToStaticMarkup(userInfo) });
                        marker.addListener('click', () => info.open(map, marker));
                    });
                });
            });
        };
    },
    componentDidMount: function() {
        this.loadMap();
        var maps_script = document.createElement("script");
        maps_script.src = "https://maps.googleapis.com/maps/api/js?key=AIzaSyD9Qqc1TdmiQ4neQy5PQoxzq-lU3PCIjqY&callback=initMap";
        maps_script.async = true;
        document.body.appendChild(maps_script);
    },


    render: function() {
        return <div id="communityMap" style={{height: "100%"}}></div>

    }
});

module.exports = Community;
