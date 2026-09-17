'use strict';

import React from 'react';
import {Container as Grid, Row, Col, Button} from 'react-bootstrap';
import {Link, withRouter} from 'react-router-dom';
import auth from './auth';
import {$c} from './Signup';
import config from './config';
import $ from 'jquery';

class ResetPassword extends React.Component {
    state = { submitted: null, success: null, error: null };
    contactFormRef = React.createRef();

    render() {
        var message;
        if (this.state.error != null) {
            message = (
				<div className="alert alert-danger">
					<p>{this.state.error}</p>
				</div>);
        } else if (this.state.success) {
            message = (
				<div className="alert alert-success">
					<p>We have sent you an email with instructions to reset your password.</p>
				</div>);
        }
        return (
            <Grid id="main-grid">
                <Row id="message">
                    {message}
                </Row>
                <Row>
                    <Col md={{ span: 4, offset: 4 }} sm={{ span: 6, offset: 3 }}>
                        <h3>Forgot your password?</h3>
                        <div>Enter your email and we&apos;ll send you a link to reset your password.</div>
                        <br/>
                        <ResetPasswordForm ref={this.contactFormRef}/>
			<Button type="button" className="btn btn-primary w-100 mt-2" onClick={this.handleSubmit}>
                            Reset Password
                        </Button>
                    </Col>
                </Row>
            </Grid>);
    }

    handleSubmit = () => {
        var showSuccess = () => {this.setState({success: true});};
        var showFailure = msg => {this.setState({error: msg});};
        var url = config.backend_url + '/accounts/password_reset/';
        if (this.contactFormRef.current && this.contactFormRef.current.isValid()) {
            var formData = this.contactFormRef.current.getFormData();
            $.post({
                url: url,
                data: formData,
                success: function (data) {
                    showFailure(data.error);
                    if (data.success) {
                        showSuccess();
                    }
                },
                error: function () {
                    showFailure('Could not complete this action');
                }
            });
        }
    };
}


class ResetPasswordForm extends React.Component {
    state = { errors: {} };
    emailRef = React.createRef();

    isValid = () => {
        var compulsoryFields = ['email'];
        var errors = {};
        compulsoryFields.forEach((field) => {
	    const node = field === 'email' ? this.emailRef.current : null;
            var value = (node && node.value ? node.value : '').trim();
            if (!value) {
                errors[field] = 'This field is required';
            }
        });
        this.setState({errors: errors});
		return Object.keys(errors).length === 0;
    };
    getFormData = () => {
        var data = {
		email: this.emailRef.current ? this.emailRef.current.value : ''
        };
        return data;
    };
    render() {
        return (
			<div className="form-horizontal">
				{this.renderTextInput('email', 'Email')}
			</div>);
    }
    renderTextInput(id, label) {
        return this.renderField(id, label,
            <input type="text" className="form-control" id={id} ref={this.emailRef}/>
        );
    }
    renderField(id, label, field) {
        return (
			<div className={$c('mb-3', {'has-error': id in this.state.errors})}>
				<label htmlFor={id} className="form-label">{label}</label>
				<div>
					{field}
				</div>
			</div>);
    }
}

class SigninInner extends React.Component {
    state = { submitted: null, success: null, successMessage: null, error: null };
    contactFormRef = React.createRef();

    render() {
        var message;
        if (this.state.error != null) {
            message = (
				<div className="alert alert-danger">
					<p>{this.state.error}</p>
				</div>);
        } else if (this.state.success && this.state.successMessage != null) {
            message = (
				<div className="alert alert-success">
					<p>{this.state.successMessage}</p>
				</div>);
        }
        return (
            <Grid id="main-grid">
                <Row id="message">
                    {message}
                </Row>
                <Row id="form">
                    <Col md={{ span: 4, offset: 4 }} sm={{ span: 6, offset: 3}}>
                        <SigninForm onSubmit={e => { this.handleSubmit(); e.preventDefault(); }} ref={this.contactFormRef}/>
                        <Button type="button" className="btn btn-primary w-100 mt-2" onClick={this.handleSubmit}>
                            Sign in
                        </Button>
                    	<div className="text-end mt-2">
                            <Link to='/reset_password'><Button variant="link">Forgot your password?</Button></Link>
                    	</div>
		    </Col>
                </Row>
            </Grid>);
    }

    handleSubmit = () => {
        if (this.contactFormRef.current && this.contactFormRef.current.isValid()) {
            var formData = this.contactFormRef.current.getFormData();
            auth.login(formData.email, formData.password, (loggedIn, error) => {
                if (loggedIn) {
                    const query = new URLSearchParams(this.props.location?.search || '');
		    var target = query.get('target');
                    if (target == null) {
                        target = '/profile';
                    }
                    this.props.history.push(target);
                } else {
		    var errorMessage = error && (error.detail || error.non_field_errors);
                    if (errorMessage === 'No active account found with the given credentials') {
                        this.setState({error: "Incorrect email/password"});
                    } else {
                        this.setState({error: errorMessage || "Could not sign in"});
                    }
                }
            });
        } else {
            this.setState({error: "Some information was missing"});
        }
    };
}

class SigninForm extends React.Component {
    state = { errors: {} };
    emailRef = React.createRef();
    passwordRef = React.createRef();

    isValid = () => {
        var compulsoryFields = ['email', 'password'];
        var errors = {};
        compulsoryFields.forEach((field) => {
	    const node = field === 'email' ? this.emailRef.current : this.passwordRef.current;
            var value = (node && node.value ? node.value : '').trim();
            if (!value) {
                errors[field] = 'This field is required';
            }
        });
        this.setState({errors: errors});
		return Object.keys(errors).length === 0;
    };
    getFormData = () => {
        var data = {
	    email: this.emailRef.current ? this.emailRef.current.value : '',
            password: this.passwordRef.current ? this.passwordRef.current.value : ''
        };
        return data;
    };
    render() {
        return (
            <form className="form" onSubmit={this.props.onSubmit}>
                {this.renderTextInput('email', 'Email')}
                {this.renderPassword('password', 'Password')}
                <input type="submit" className="d-none" />
            </form>);
    }
    renderTextInput(id, label) {
        return this.renderField(id, label,
            <input type="text" className="form-control" id={id} ref={this.emailRef}/>
        );
    }
    renderPassword(id, label) {
        return this.renderField(id, label,
            <input type="password" className="form-control" id={id} ref={this.passwordRef}/>
        );
    }
    renderField(id, label, field) {
        return (
			<div className={$c('mb-3', {'has-error': id in this.state.errors})}>
				<label htmlFor={id} className="form-label">{label}</label>
				<div>
				    {field}
				</div>
			</div>);
    }
}

const Signin = withRouter(SigninInner);
export { Signin, ResetPassword };
export default { Signin, ResetPassword };
