'use strict';

var React = require('react');
var PureRenderMixin = require('./PureRenderMixin'); // deep-equals version of PRM
var {Grid, Col, Row, Button, Table} = require('react-bootstrap');
var backend = require('backend');
var config  = require('./config');
var {Role} = require('./Signup');
var {Navigation, Link} = require('react-router');
var {Pagination} = require('react-data-components-bd2k');
var _ = require('underscore');
var placeholder = require('./img/placeholder.png');
var auth = require('./auth');
var Rx = require('rx');

var Community = React.createClass({
    mixins: [PureRenderMixin, Navigation],
    componentWillMount: function () {
        this.fetch(this.state);
    },
    getInitialState: function () {
        return {
            data: [],
            search: "",
            pageLength: 10,
            page: 0,
            totalPages: 1
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
    onChangeSearch: function (search) {
        this.setStateFetch({search: search, page: 0});
    },
    onFilterRole: function (role) {
        if (Role.get(role).length === 3) {
            this.refs['community-search'].appendSearch(Role.get(role)[2]);
        }
    },
    logout: function() {
        auth.logout();
        this.forceUpdate();
    },
    componentDidMount() {
        var searchq = this.searchq = new Rx.Subject();
        this.subs = searchq.debounce(500).subscribe(this.onChangeSearch);
    },
    componentWillUnmount() {
        this.subs.dispose();
    },
    render: function () {
        var queryParams = this.context.router.getCurrentQuery();
        var message;
        if (queryParams.registrationSuccess === "true") {
            message = (
				<div className="alert alert-success">
					<p>Thanks for signing up. We have sent you an email with a confirmation link to complete your registration.</p>
				</div>);
        } else if (queryParams.updateSuccess === "true") {
            message = (
                <div className="alert alert-success">
                    {queryParams.subscribe === "true" &&
                        <p>Your profile has been edited successfully, and you have been added to the mailing list.</p>}
                    {queryParams.subscribe === "false" &&
                        <p>Your profile has been edited successfully, and you have been removed from the mailing list.</p>}
                    {queryParams.subscribe === undefined &&
                        <p>Your profile has been edited successfully.</p>}
                </div>);
        }
        var {data, page, totalPages, error} = this.state;
        var rows = _.map(data, row => {

            var avatar;
            if (row.has_image) {
                var avatarLink = config.backend_url + '/site_media/media/' + row.id;
                avatar = <img className="avatar" src={avatarLink} type="image/image"/>;
            } else {
                avatar = <img src={placeholder}/>;
            }

            var {city, state, country} = row;
            var locationString = _.values(_.pick({city, state, country}, v => v)).join(', ');

            return (<tr>
                <td width="120px">
                    {avatar}
                </td>
                <td>
                    <span id="name"><h3>{row.firstName} {row.lastName}{row.title.length ? ", " + row.title : ""}</h3></span>
                    <span id="role"><h4>{Role.other(row.role) ? row["role_other"] : Role.get(row.role)[2]}</h4></span>
                    <span id="institution">{row.institution}</span>
                    <span id="location">{locationString}</span>
                    <span id="contact">{row.email} {row["phone_number"]}</span>
                </td>
            </tr>);
        });

        return (error ? <p>{error}</p> :
            <Grid id="main-grid">
                <Row id="message"> {message} </Row>
                <Row>
                    <Col sm={10} smOffset={1} md={8} mdOffset={2}>
                        <span>The BRCA Exchange supports the exchange of information about BRCA1 and BRCA2 variants. Show your support by joining our mailing list and/or listing your name below as one of our supporters.</span>
                    </Col>
                </Row>
                <Row>
                    <Col md={8} mdOffset={2} sm={10} smOffset={1}>
                        <CommunityMap onFilterRole={this.onFilterRole} search={this.state.search}/>
                    </Col>
                </Row>
                <Row>
                    <Col className="text-center" sm={12} smOffset={0} md={10} mdOffset={1}>
                        <Link to="/mailinglist"><Button disabled={config.environment === 'beta' && 'disabled'}>Join the mailing list only</Button></Link>&nbsp;
                        <Link to="/signup"><Button disabled={config.environment === 'beta' && 'disabled'}>Join the mailing list and add me to the supporters</Button></Link>&nbsp;
                        <Link to="/signin"><Button disabled={config.environment === 'beta' && 'disabled'}>Edit your profile</Button></Link>
                    </Col>
                </Row>
                <Row>
                    <Col className="text-center">
                        <h2>We Support the BRCA Exchange</h2>
                    </Col>

                </Row>
                <Row>
                    <Col className="btm-buffer" sm={10} smOffset={1} md={8} mdOffset={2}>
                        <Col sm={6} lg={5} style={{paddingRight: "0"}}>
                            <h4>Search for a community member:</h4>
                        </Col>
                        <Col sm={6} lg={7}>
                            <CommunitySearch ref="community-search" onChange={s => this.searchq.onNext(s)}/>
                        </Col>
                    </Col>
                </Row>
                <Row>
                    <Col md={8} mdOffset={2} sm={10} smOffset={1}>
                        <Table className="community" striped bordered>
                            <tbody>
                                {rows}
                            </tbody>
                        </Table>
                    </Col>
                </Row>
                <Row>
                    <Col sm={10} smOffset={1} md={8} mdOffset={2}>
                        <Pagination
                            className="pagination pull-right-sm"
                            currentPage={page}
                            totalPages={totalPages}
                            onChangePage={this.onChangePage} />
                    </Col>
                </Row>
            </Grid>
        );
    }
});

var CommunityMap = React.createClass({
    getInitialState: () => ({ roles: [] }),
    shouldComponentUpdate({search}) {
        return this.props.search !== search;
    },
    onFilterRole: function(role) {
        var roles = this.state.roles.slice();
        roles[role] = !roles[role];
        this.setState({ roles: roles });
        this.filterSub.onNext(role);
    },
    componentDidMount: function() {

        var initMap = function () {
            var self = this;
            var infowindow;
            var markers = this.markers = [];
            var map = this.map = new google.maps.Map(document.getElementById('communityMap'), {
                center: {lat: 17, lng: -2.5},
                zoom: 1,
                scrollwheel: false,
                mapTypeControl: false,
                streetViewControl: false,
                styles: [{
                    "featureType": "administrative",
                    "elementType": "geometry.fill",
                    "stylers": [{ "visibility": "off" }]
                }]
            });

            map.addListener('click', function() {
                if (infowindow) {
                    infowindow.close();
                }
                infowindow = undefined;
            });

            var legendControl = document.createElement('div');
            var roleMarkers = Role.options.map(role =>
                <div className="map-legend-col" data-role={role[0]}>
                    <img src={require(`./img/map/${role[0]}key.png`)} /> {role.length === 3 ? role[2] : role[1]}
                </div>
            );
            roleMarkers = _.values(_.groupBy(roleMarkers, (item, index) => index % 4)).map(group => <div className="map-legend-row">{group}</div>);
            legendControl.innerHTML = React.renderToStaticMarkup(
                <div>
                    <img src={require("./img/map/1.png")} className="map-legend-icon" />
                    <div className="map-legend-slide">
                        <div className="map-legend-full">
                            {roleMarkers}
                        </div>
                    </div>
                </div>
            );

            Array.prototype.slice.call(legendControl.querySelectorAll(".map-legend-col")).map(function(el) {
                el.addEventListener('click', function() {
                    self.onFilterRole(this.attributes['data-role'].value);
                    this.classList.toggle('selected');
                });
            });
            legendControl.className = "map-legend";
            legendControl.index = 1;
            map.controls[google.maps.ControlPosition.TOP_RIGHT].push(legendControl);

            this.updateMap = (function() { //eslint-disable-line no-extra-parens
                var roleFilter = [...this.state.roles.keys()].filter(index => this.state.roles[index]);
                backend.userLocations(this.props.search, roleFilter).subscribe(({data}) => {
                    var newMarkers = [];
                    markers.map(m => {
                        var idx = data.findIndex(({id}) => m.userID === id);
                        if (idx !== -1) {
                            newMarkers.push(m);
                            data.splice(idx, 1);
                        } else {
                            m.setMap(null);
                        }
                    });
                    markers = self.markers = newMarkers;
                    /*eslint-disable camelcase*/
                    _.map(data, ({id, firstName, lastName, title, role, role_other, institution, latitude, longitude, has_image})  => {
                        if (latitude !== "" || longitude !== "") {
                            var avatar;
                            if (has_image) {
                                let avatar_link = config.backend_url + '/site_media/media/' + id;
                                avatar = <object className="avatar" data={avatar_link} type="image/jpg"/>;
                            } else {
                                avatar = <img className="avatar" src={placeholder}/>;
                            }
                            var userInfo = (<div className="map-info-window">
                                {avatar}
                                <div>
                                    <span>{firstName} {lastName}{title.length ? "," : ""} {title}</span><br />
                                    <span id="role">{Role.other(role) ? role_other : Role.get(role)[2]}</span><br />
                                    <span>{institution}</span>
                                </div>
                            </div>);
                            /*eslint-enable camelcase*/
                            var marker = new google.maps.Marker({
                                position: { lat: parseFloat(latitude), lng: parseFloat(longitude) },
                                map: map,
                                title: `${firstName} ${lastName}${title.length ? "," : ""} ${title}`,
                                icon: {
                                    url: require(`./img/map/${role}.png`),
                                    scaledSize: new google.maps.Size(20, 32)
                                }
                            });
                            marker.userID = id;
                            markers.push(marker);
                            var info = new google.maps.InfoWindow({content: React.renderToStaticMarkup(userInfo) });
                            marker.addListener('click', () => {
                                if (infowindow) {
                                    infowindow.close();
                                }
                                infowindow = info;
                                info.open(map, marker);
                            });
                        }
                    });
                });
                if (infowindow) {
                    infowindow.close();
                }
                infowindow = undefined;
            }).bind(this);
            this.updateMap();
            var filterSub = this.filterSub = new Rx.Subject();
            this.subs = filterSub.debounce(500).subscribe(this.updateMap);
        };


        google.load('maps', '3', {callback: initMap.bind(this), "other_params": "key=" + config.maps_key});
    },

    componentDidUpdate: function() {
        this.updateMap();
    },
    componentWillUnmount: function() { window.removeEventListener('resize', this.handleResize); },

    render: function() {
        return <div id="communityMap"></div>;
    }
});

var CommunitySearch = React.createClass({
    getInitialState: () => ({ search: "" }),
    onChange: function(e) {
        this.setState({search: e.target.value});
        this.props.onChange(e.target.value);
    },
    onSubmit: function (ev) {
        ev.preventDefault();
    },
    appendSearch: function(term) {
        this.setState({search: `${this.state.search.trim()} ${term}`.trim()});
        this.props.onChange(this.state.search);
    },
    render: function() {
        return (<div className='search-box'>
            <form onSubmit={this.onSubmit} style={{display: 'inline'}}>
                <input type='submit' className='input-sm'style={{display: 'none'}} />
                <div className='text-nowrap help-target'>
                    <div>
                        <input className='community-search-input' placeholder='name, organization, city, etc.' type='text'
                            onChange={this.onChange} value={this.state.search} />
                    </div>
                </div>
            </form>
        </div>);
    }
});

module.exports = Community;
