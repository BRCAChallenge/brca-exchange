
import React from 'react';
import { Grid, Row } from 'react-bootstrap';
import { State, Navigation } from 'react-router';
import $ from 'jquery';
import config from '../config';


var ConfirmEmail = React.createClass({
    mixins: [State, Navigation],
    getInitialState: function () {
        return {
            success: null
        };
    },
	activate: function(result) {
		this.setState({success: result.success});
	},
    componentDidMount: function () {
        var activationCode = this.context.router.getCurrentParams().activationCode;
        var url = config.backend_url + '/accounts/confirm/' + activationCode + '/';
        $.get(url, this.activate);
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


export default ConfirmEmail;
