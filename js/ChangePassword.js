'use strict';

var React = require('react');
var {Grid, Row, Col, Button} = require('react-bootstrap');
var {State, Navigation} = require('react-router');
var {trim, $c} = require('./Signup');

var $ = require('jquery');
var config = require('./config');


var ChangePassword = React.createClass({
    mixins: [State, Navigation],
    getInitialState: function () {
        return {
            success: null
        }
    },
    componentDidMount: function () {
        var showFailure = (msg => {this.setState({error: msg})});
        var resetToken = this.context.router.getCurrentParams().resetToken;
        var url = config.backend_url + '/accounts/update_password/' + resetToken + '/';
        $.post({
            url: url,
            success: function (data) {
                if (data.invalid_token) {
                    this.setState({'invalid_token': true})
                }
            }.bind(this)
        });
        },
    render: function () {
        // If the token is invalid, show an error and don't show the form.
        if (this.state.invalid_token) {
            return <Grid id="main-grid"> <Row id="message">
                <div className="alert alert-danger">
                    <p>Invalid link</p>
                </div>
            </Row> </Grid>
        }

        var message;
       if (this.state.success) {
            message = <div className="alert alert-success">
                <p>Your password has been updated. You can now sign in using it.</p>
            </div>
        }
        return (
            <Grid id="main-grid"> <Row id="message">
                {message}
            </Row>
                <Row>
                    <Col md={8} mdOffset={3}>
                        <h3>Create a new password</h3>
                    </Col>
                </Row>
                <Row id="form">
                    <Col md={8} mdOffset={2}>
                        <ChangePasswordForm ref="contactForm"/>
                    </Col>
                </Row>
                <Row id="submit">
                    <Col md={6} mdOffset={3}>
                        <Button type="button" className="btn btn-primary btn-block" onClick={this.handleSubmit}>
                            Save
                        </Button>
                    </Col>
                </Row>
            </Grid>);
    },
    handleSubmit: function () {
        var showSuccess = (() => {this.setState({success: true})});
        var showFailure = (msg => {this.setState({error: msg})});
        var resetToken = this.context.router.getCurrentParams().resetToken;

        var url = config.backend_url + '/accounts/update_password/' + resetToken + '/';
        if (this.refs.contactForm.isValid()) {
            var formData = this.refs.contactForm.getFormData();
            $.post({
                url: url,
                data: formData,
                success: function (data) {
                    showFailure(data.error);
                    if (data.success) {
                        showSuccess()
                    }
                },
                error: function (xhr, ajaxOptions, thrownError) {
                    showFailure('Could not complete this action');
                }.bind(this)
            });
        }
    }
});


var ChangePasswordForm = React.createClass({
    getInitialState: function () {
        return {errors: {}}
    },
    isValid: function () {
        var compulsory_fields = ['password', 'password_confirm'];
        var errors = {};
        if (this.refs.password.getDOMNode().value != this.refs.password_confirm.getDOMNode().value) {
            errors["password_confirm"] = "The passwords don't match"
        }
        compulsory_fields.forEach(function (field) {
            var value = trim(this.refs[field].getDOMNode().value)
            if (!value) {
                errors[field] = 'This field is required'
            }
        }.bind(this));
        this.setState({errors: errors});
        var isValid = true;
        for (var error in errors) {
            isValid = false;
            break;
        }

        return isValid
    },
    getFormData: function () {
        var data = {
            password: this.refs.password.getDOMNode().value
            , password_confirm: this.refs.password_confirm.getDOMNode().value
        };
        return data
    },
    render: function () {
        return <div className="form-horizontal">
            {this.renderPassword('password', 'Password')}
            {this.renderPassword('password_confirm', 'Confirm Password')}
        </div>
    },
    renderPassword: function (id, label) {
        return this.renderField(id, label,
            <input type="password" className="form-control" id={id} ref={id}/>
        )
    },
    renderField: function (id, label, field) {
        return <div className={$c('form-group', {'has-error': id in this.state.errors})}>
            <label htmlFor={id} className="col-sm-4 control-label">{label}</label>
            <div className="col-sm-6">
                {field}
            </div>
        </div>
    }
});


module.exports = ({
    ChangePassword: ChangePassword
});
