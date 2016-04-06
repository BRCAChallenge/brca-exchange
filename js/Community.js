'use strict';

var React = require('react');
var PureRenderMixin = require('./PureRenderMixin'); // deep-equals version of PRM
var {Grid, Col, Row, Button, Table} = require('react-bootstrap');
var backend = require('backend');
var {Navigation, Link} = require('react-router');
var {Pagination} = require('react-data-components-bd2k');
var _ = require('underscore');
var placeholder = require('./img/placeholder.png');

var Community = React.createClass({
    mixins: [PureRenderMixin, Navigation],
    componentWillMount: function () {
        backend.users(this.state).subscribe(
            resp => this.setState({...resp}), // set data, count, totalPages
            () => this.setState({error: 'Problem connecting to server'}));
    },
    getInitialState: function () {
        return {
            data: [],
            pageLength: 10,
            page: 0,
            totalPages: 0,
            windowWidth: window.innerWidth
        };
    },
    fetch: function (state) {
        var {pageLength, page} = state;
        this.fetchq.onNext({pageLength,page});
    },
    setStateFetch: function (opts) {
        var newState = {...this.state, ...opts};
        this.setState(newState);
        this.fetch(newState);
    },
    onChangePage: function (pageNumber) {
        this.setStateFetch({page: pageNumber});
    },
    render: function () {
        var {data, page, totalPages, error} = this.state;

        var rows = _.map(data, row => {

            var avatar;
            if (row.has_image) {
                var avatar_link = backend.databaseUrl + '/site_media/media/' + row['id']
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
                <Row>
                    <div className='text-center Variant-detail-title'>
                        <h3>BRCA Community</h3>
                    </div>
                </Row>
                <Row>
                    <Col sm={4}>
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
                    <Col md={8} mdOffset={2}>
                        <Link to="/signup"><Button>Join our mailing list and this community space</Button></Link>
                    </Col>
                </Row>
            </Grid>
        );
    }
});

module.exports = Community;
