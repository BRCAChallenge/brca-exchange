'use strict';

var React = require('react');
var {Grid, Row} = require('react-bootstrap');
var {State, Navigation} = require('react-router');
var $ = require('jquery');
var config = require('./config');


var ConfirmEmail = React.createClass({
    mixins: [State, Navigation],
    getInitialState: function () {
        return {
            success: null
        }
    },
    componentDidMount: function () {
        var activationCode = this.context.router.getCurrentParams().activationCode;
        var url = config.backend_url + '/accounts/confirm/' + activationCode + '/';
        $.get(url, function (result) {
            this.setState({success: result.success});
        }.bind(this));
    },
    render: function () {
        return (
            <Grid id="main-grid">
                {this.state.success &&
                <Row id="message" className="alert alert-success">
                    Thanks for confirming your email.
                </Row>
                }

                {!this.state.success &&
                <Row id="message" className="alert alert-danger">
                    Sorry, this activation link is not valid. </Row>
                }
            </Grid>);
    }
});


module.exports = ({
    ConfirmEmail: ConfirmEmail
});
