'use strict';

import React from 'react';
import ReactDOMServer from 'react-dom/server';
import { Container as Grid, Col, Row, Button, Table } from 'react-bootstrap';
import backend from './backend';
import config from './config';
import { Role } from './Signup';
import { Link, withRouter } from 'react-router-dom';
// TODO: Uncomment when react-data-components-brcaex is updated/replaced
// import { Pagination } from 'react-data-components-brcaex';
import _ from 'underscore';
import placeholder from './img/placeholder.png';
import auth from './auth';

// RxJS imports
import { Subject } from 'rxjs';
import { debounceTime } from 'rxjs/operators';

// ---- Google Maps loader (async, no jsapi) ----
function loadGoogleMaps(key, cb) {
  if (!key) {
    console.warn('Maps key missing; skipping map init');
    return;
  }
  if (window.google && window.google.maps) {
    cb();
    return;
  }
  const script = document.createElement('script');
  // async load + callback to your init function
  window.__initMap = cb;
  script.src = `https://maps.googleapis.com/maps/api/js?key=${encodeURIComponent(
    key
  )}&callback=__initMap&loading=async`;
  script.async = true;
  script.defer = true;
  document.head.appendChild(script);
}

class Community extends React.Component {
    constructor(props) {
	super(props);
	this.state = {
            data: [],
            search: "",
            pageLength: 10,
            page: 0,
            totalPages: 1
        };
	this.communitySearchRef = React.createRef();
    }

    setPages = ({data, count}) => {
        return {
            data,
            count,
            totalPages: Math.ceil(count / this.state.pageLength)
        };
    };

    fetch = (state) => {
        backend.users(state).subscribe(
            resp => this.setState(this.setPages(resp)), // set data, count, totalPages
            () => this.setState({error: 'Problem connecting to server'}));
    };

    setStateFetch = (opts) => {
        var newState = {...this.state, ...opts};
        this.setState(newState);
        this.fetch(newState);
    };

    onChangePage = (pageNumber) => {
        this.setStateFetch({page: pageNumber});
    };

    onChangeSearch = (search) => {
        this.setStateFetch({search: search, page: 0});
    };

    onFilterRole = (role) => {
        if (Role.get(role).length === 3 && this.communitySearchRef.current) {
            this.communitySearchRef.current.appendSearch(Role.get(role)[2]);
        }
    };

    logout = () => {
        auth.logout();
        this.forceUpdate();
    };

    componentDidMount() {
	this.fetch(this.state);
        var searchq = this.searchq = new Subject();
        this.subs = searchq.pipe(
            debounceTime(500)
        ).subscribe(this.onChangeSearch);
    };

    componentWillUnmount() {
        this.subs.unsubscribe();
    }

    render() {
        const query = new URLSearchParams(this.props.location?.search || "");
        var message;
        if (query.get("registrationSuccess") === "true") {
            message = (
				<div className="alert alert-success">
					<p>Thanks for signing up. We have sent you an email with a confirmation link to complete your registration. After you complete your registration our administrator will confirm your profile and it will appear on the Community Pages.</p>
				</div>);
        } else if (query.get("updateSuccess") === "true") {
            message = (
                <div className="alert alert-success">
                    {query.get("subscribe") === "true" &&
                        <p>Your profile has been edited successfully, and you have been added to the mailing list.</p>}
                    {query.get("subscribe") === "false" &&
                        <p>Your profile has been edited successfully, and you have been removed from the mailing list.</p>}
                    {!query.has("subscribe") &&
                        <p>Your profile has been edited successfully.</p>}
                </div>);
        }
        var {data, page, totalPages, error} = this.state;
        var rows = _.map(data, (row, i) => {

            var avatar;
            if (row.has_image) {
                var avatarLink = config.backend_url + '/site_media/media/' + row.id;
                avatar = <img className="avatar" src={avatarLink} type="image/image"/>;
            } else {
                avatar = <img src={placeholder}/>;
            }

            var {city, state, country} = row;
            var locationString = _.values(_.pick({city, state, country}, v => v)).join(', ');

            return (<tr key={row.id ?? i}>
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
                    <Col md={{ span: 10, offset: 1 }} sm={12}>
                        <p className="community-message">The BRCA Exchange supports the exchange of information about BRCA1 and BRCA2 variants. Show your support by joining our global community!  By showing your support, you will help us demonstrate the value of this resource, which will help keep it freely available to all.</p>
                    </Col>
                </Row>
                <Row>
                    <Col md={{ span: 10, offset: 1 }} sm={12}>
                        <CommunityMap onFilterRole={this.onFilterRole} search={this.state.search}/>
                    </Col>
                </Row>
                <Row>
                    <Col className="text-center" md={{ span: 10, offset: 1 }} sm={12}>
                        <Link to="/signup"><Button disabled={config.environment === 'beta'}>Join the community</Button></Link>&nbsp;
                        <p className="community-disclaimer">We will add you to the BRCA Exchange News mailing list. You can unsubscribe at any time.</p>
                        <p>To update or remove your profile, please <a href="mailto:brca-exchange-contact@genomicsandhealth.org?subject=Update Personal Information">contact us</a>.</p>
                    </Col>
                </Row>
                <Row>
                    <Col className="text-center">
                        <h2>We Support the BRCA Exchange</h2>
                    </Col>

                </Row>
                <Row>
                    <Col className="btm-buffer" md={{ span: 10, offset: 1 }} sm={12}>
                        <Row>
			    <Col sm={6} lg={5} style={{paddingRight: "0"}}>
                        	<h4>Search for a community member:</h4>
                            </Col>
                            <Col sm={6} lg={7}>
                        	<CommunitySearch ref={this.communitySearchRef} onChange={s => this.searchq.next(s)}/>
                            </Col>
			</Row>
                    </Col>
                </Row>
                <Row>
                    <Col md={{ span: 10, offset: 1 }} sm={12}>
                        <Table className="community" striped bordered>
                            <tbody>
                                {rows}
                            </tbody>
                        </Table>
                    </Col>
                </Row>
                <Row>
                    <Col md={{ span: 10, offset: 1 }} sm={12}>
                        <p style={{verticalAlign: 'bottom', display: 'inline-block'}}>
                        {`${(this.state.page * this.state.pageLength) + 1}-${Math.min((this.state.page + 1) * this.state.pageLength, this.state.count)} out of ${this.state.count} members`}
                        </p>

                        {/* TODO: Uncomment when react-data-components-brcaex is updated/replaced */}
                        {/*
                        <Pagination
                            className="pagination pull-right-sm"
                            currentPage={page}
                            totalPages={totalPages}
                            onChangePage={this.onChangePage} />
                        */}

                        {/* TEMPORARY: Basic pagination buttons until Pagination component is available */}
                        <div className="pagination pull-right-sm" style={{display: 'inline-block', marginLeft: '20px'}}>
                            <Button
                                variant="secondary"
                                size="sm"
                                disabled={page === 0}
                                onClick={() => this.onChangePage(page - 1)}
                                style={{marginRight: '5px'}}
                            >
                                Previous
                            </Button>
                            <span style={{padding: '0 10px'}}>
                                Page {page + 1} of {totalPages}
                            </span>
                            <Button
                                variant="secondary"
                                size="sm"
                                disabled={page >= totalPages - 1}
                                onClick={() => this.onChangePage(page + 1)}
                            >
                                Next
                            </Button>
                        </div>
                    </Col>
                </Row>
            </Grid>
        );
    }
}

class CommunityMap extends React.Component {
    constructor(props) {
	super(props);
	this.state = { roles: [] };
    }

    shouldComponentUpdate({search}) {
        return this.props.search !== search;
    }

    onFilterRole = (role) => {
        var roles = this.state.roles.slice();
        roles[role] = !roles[role];
        this.setState({ roles: roles });
        this.filterSub.next(role);
    };

    componentDidMount() {

        var initMap = function () {
            var self = this;
            var infowindow;
            var markers = this.markers = [];
            var map = this.map = new google.maps.Map(document.getElementById('communityMap'), {
                center: {lat: 17, lng: 9.4},
                zoom: 2,
                scrollwheel: false,
                mapTypeControl: true,
                mapTypeId: 'roadmap',
                streetViewControl: false,
                styles: [{
                    "featureType": "administrative",
                    "elementType": "geometry.fill",
                    "stylers": [{ "visibility": "off" }]
                }],
                zoomControlOptions: {
                    position: google.maps.ControlPosition.LEFT_BOTTOM
                }
            });

            map.addListener('click', function() {
                if (infowindow) {
                    infowindow.close();
                }
                infowindow = undefined;
            });

            var legendControl = document.createElement('div');
            var roleMarkers = Role.options.map(role =>
                (<div key={`role-${role[0]}`} className="map-legend-col" data-role={role[0]}>
                    <img src={require(`./img/map/${role[0]}key.png`)} /> {role.length === 3 ? role[2] : role[1]}
                </div>)
            );
            roleMarkers = _.values(_.groupBy(roleMarkers, (item, index) => index % 4)).map((group, i) => <div key={`legend-row-${i}`} className="map-legend-row">{group}</div>);
	    legendControl.innerHTML = ReactDOMServer.renderToStaticMarkup(
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
                                    scaledSize: new google.maps.Size(15, 24)
                                }
                            });
                            marker.userID = id;
                            markers.push(marker);
                            var info = new google.maps.InfoWindow({content: ReactDOMServer.renderToStaticMarkup(userInfo) });
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
            var filterSub = this.filterSub = new Subject();
            this.subs = filterSub.pipe(
                debounceTime(500)
            ).subscribe(this.updateMap);
        };

	loadGoogleMaps(config.maps_key, initMap.bind(this));
    }

    componentDidUpdate() {
        this.updateMap();
    }

    componentWillUnmount() {
        window.removeEventListener('resize', this.handleResize);
        if (this.subs) {
            this.subs.unsubscribe();
        }
    }

    render() {
        return <div id="communityMap" />;
    }
}

class CommunitySearch extends React.Component {
    constructor(props) {
	super(props);
	this.state = {
	    search: "",
            placeholder: "name, organization, city, etc."
	};
    }

    onChange = (e) => {
        this.setState({search: e.target.value});
        this.props.onChange(e.target.value);
    };

    onSubmit = (ev) => {
        ev.preventDefault();
    };

    onBlur = () => {
        this.setState({placeholder: "name, organization, city, etc."});
    };

    onFocus = () => {
        this.setState({placeholder: ""});
    };

    appendSearch = (term) => {
        this.setState({search: `${this.state.search.trim()} ${term}`.trim()});
        this.props.onChange(this.state.search);
    };

    render() {
        return (<div className='search-box'>
            <form onSubmit={this.onSubmit} style={{display: 'inline'}}>
                <input type='submit' className='input-sm'style={{display: 'none'}} />
                <div className='text-nowrap help-target'>
                    <div>
                        <input className='community-search-input'
                               onBlur={this.onBlur}
                               onFocus={this.onFocus}
                               placeholder={this.state.placeholder}
                               type='text'
                               onChange={this.onChange}
                               value={this.state.search} />
                        <span className="fa fa-search search-box-icon"/>
                    </div>
                </div>
            </form>
        </div>);
    }
}

export default withRouter(Community);
