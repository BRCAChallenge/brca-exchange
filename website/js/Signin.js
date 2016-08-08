'use strict';

var React = require('react');
var content = require('./content');
var RawHTML = require('./RawHTML');
var {Grid, Row, Col, Button} = require('react-bootstrap');
var {State, Navigation, Link} = require('react-router');
var auth = require('./auth');
var {trim, $c} = require('./Signup');
var config  = require('./config')
var $ = require('jquery');

var ResetPassword = React.createClass({
    mixins: [State, Navigation],
    getInitialState: function () {
        return {
            submitted: null,
            success: null
        }
    },
    render: function () {
        var message;
        if (this.state.error != null) {
            message = <div className="alert alert-danger">
                <p>{this.state.error}</p>
            </div>
        } else if (this.state.success) {
            message = <div className="alert alert-success">
                <p>We have sent you an email with instructions to reset your password.</p>
            </div>
        }
        return (
            <Grid id="main-grid">
                <Row id="message">
                    {message}
                </Row>
                <Row>
                    <Col md={8} mdOffset={3}>
                        <h3>Forgot your password?</h3>
                        <div>Enter your email and we'll send you a link to reset your password.</div>
                        <br/>
                    </Col>
                </Row>
                <Row id="form">
                    <Col md={8} mdOffset={2}>
                        <ResetPasswordForm ref="contactForm"/>
                    </Col>
                </Row>
                <Row id="submit">
                    <Col md={6} mdOffset={3}>
                        <Button type="button" className="btn btn-primary btn-block" onClick={this.handleSubmit}>
                            Reset Password
                        </Button>
                    </Col>
                </Row>
            </Grid>)
    },

    handleSubmit: function () {
        var showSuccess = (() => {this.setState({success: true})});
        var showFailure = (msg => {this.setState({error: msg})});
        var url = config.backend_url + '/accounts/password_reset/';
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


var ResetPasswordForm = React.createClass({
    getInitialState: function () {
        return {errors: {}}
    },
    isValid: function () {
        var compulsory_fields = ['email'];
        var errors = {};
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
            email: this.refs.email.getDOMNode().value
        };
        return data
    },
    render: function () {
        return <div className="form-horizontal">
            {this.renderTextInput('email', 'Email')}
        </div>
    },
    renderTextInput: function (id, label) {
        return this.renderField(id, label,
            <input type="text" className="form-control" id={id} ref={id}/>
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

var Signin = React.createClass({
    mixins: [State, Navigation],
    getInitialState: function () {
        return {
            submitted: null,
            success: null,
            successMessage: null
        }
    },
    render: function () {
        var message;
        if (this.state.error != null) {
            message = <div className="alert alert-danger">
                <p>{this.state.error}</p>
            </div>
        } else if (this.state.success && this.state.successMessage != null) {
            message = <div className="alert alert-success">
                <p>{this.state.successMessage}</p>
            </div>
        }
        return (
            <Grid id="main-grid">
                <Row id="message">
                    {message}
                </Row>
                <Row id="form">
                    <Col md={8} mdOffset={2}>
                        <SigninForm onSubmit={e => { this.handleSubmit(); e.preventDefault() }} ref="contactForm"/>
                    </Col>
                </Row>
                <Row id="submit">
                    <Col md={8} mdOffset={2}>
                            <div className="form-group" style={{marginLeft: "-15px", marginRight: "-15px"}}> 
                                <label className="col-sm-4 control-label"></label>
                                <Col sm={6}>
                                    <Button type="button" className="btn btn-primary btn-block" onClick={this.handleSubmit}>
                                        Sign in
                                    </Button>
                                </Col>
                            </div>
                    </Col>
                </Row>
                <Row id="submit">
                    <Col sm={10} md={6} mdOffset={3}>
                        <Link className="pull-right" to='/reset_password'><Button bsStyle="link">Forgot your password?</Button></Link>
                    </Col>
                </Row>
            </Grid>)
    },

    handleSubmit: function () {
        if (this.refs.contactForm.isValid()) {
            var formData = this.refs.contactForm.getFormData();
            auth.login(formData.email, formData.password, (loggedIn, error) => {
                if (loggedIn) {
                    var target = this.getQuery()["target"];
                    if (target == null) {
                        target = '/profile';
                    }
                    this.transitionTo(target)
                } else {
                    if (error.non_field_errors == 'User account is disabled.') {
                        var showSuccess = (() => {this.setState({success: true, successMessage: "Activation email sent."})});
                        var showFailure = (msg => {this.setState({error: msg})});
                        var resend_activation = function() {
                            $.post({
                                url: config.backend_url + "/accounts/resend-activation/",
                                data: {email: formData.email},
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
                        };
                        var activation_message =
                            <span>
                                This account has not yet been activated. Please check your email for an activation link, or <a href="#" onClick={resend_activation}>resend activation</a>.
                            </span>
                        this.setState({error: activation_message});
                    } else if (error.non_field_errors == 'Unable to login with provided credentials.') {
                        this.setState({error: "Incorrect email/password"});
                    } else {
                        this.setState({error: error.non_field_errors});
                    }
                }
            });
        } else {
            this.setState({error: "Some information was missing"});
        }
    }
});

var SigninForm = React.createClass({
    getInitialState: function () {
        return {errors: {}}
    },

    isValid: function () {
        var compulsory_fields = ['email', 'password'];
        var errors = {};
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
            email: this.refs.email.getDOMNode().value
            , password: this.refs.password.getDOMNode().value
        };
        return data
    },
    render: function () {
        return <form className="form-horizontal" onSubmit={this.props.onSubmit}>
            {this.renderTextInput('email', 'Email')}
            {this.renderPassword('password', 'Password')}
            <input type="submit" className="hidden" />
        </form>
    },
    renderTextInput: function (id, label) {
        return this.renderField(id, label,
            <input type="text" className="form-control" id={id} ref={id}/>
        )
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
    Signin: Signin,
    ResetPassword: ResetPassword
});
