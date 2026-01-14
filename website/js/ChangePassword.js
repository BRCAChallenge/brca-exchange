'use strict';

import React from 'react';
import ReactDOM from 'react-dom';
import { Container as Grid, Row, Col, Button } from 'react-bootstrap';
import { withRouter } from 'react-router-dom';
import { $c } from './Signup';
import $ from 'jquery';
import config from './config';

class ChangePasswordInner extends React.Component {
    state = { success: null, invalid_token: false, error: null };
    contactFormRef = React.createRef();

    receiveToken = (data) => {
        if (data && data.invalid_token) {
            this.setState({ invalid_token: true });
        }
    };

    componentDidMount() {
        const resetToken =
            (this.props.match && this.props.match.params && this.props.match.params.resetToken) ||
            (this.props.params && this.props.params.resetToken);
        var url = config.backend_url + '/accounts/check_password_token/' + resetToken + '/';
        $.post({
            url: url,
            success: this.receiveToken
        });
    }

    render() {
        // If the token is invalid, show an error and don't show the form.
        if (this.state.invalid_token) {
            return (
				<Grid id="main-grid"> <Row id="message">
					<div className="alert alert-danger">
						<p>Invalid link</p>
					</div>
				</Row> </Grid>);
        }

       var message;
       if (this.state.error != null) {
            message = (
                <div className="alert alert-danger">
                    <p>{this.state.error}</p>
                </div>
            );
	} else if (this.state.success) {
            message = (
				<div className="alert alert-success">
					<p>Your password has been updated. You can now sign in using it.</p>
				</div>);
        }
        return (
            <Grid id="main-grid"> <Row id="message">
                {message}
            </Row>
                <Row>
                    <Col md={{ span: 8, offset: 3 }}>
                        <h3>Create a new password</h3>
                    </Col>
                </Row>
                <Row id="form">
                    <Col md={{ span: 8, offset: 2 }}>
                        <ChangePasswordForm ref={this.contactFormRef}/>
                    </Col>
                </Row>
                <Row id="submit">
                    <Col md={{ span: 6, offset: 3 }}>
                        <Button type="button" className="btn btn-primary btn-block" onClick={this.handleSubmit}>
                            Save
                        </Button>
                    </Col>
                </Row>
            </Grid>);
    }

    handleSubmit = () => {
        var showSuccess = () => {this.setState({success: true});};
        var showFailure = msg => {this.setState({error: msg});};
        const resetToken =
            (this.props.match && this.props.match.params && this.props.match.params.resetToken) ||
            (this.props.params && this.props.params.resetToken);

        var url = config.backend_url + '/accounts/update_password/' + resetToken + '/';
	const form = this.contactFormRef.current;
        if (form && form.isValid()) {
            var formData = form.getFormData();
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


class ChangePasswordForm extends React.Component {
    state = { errors: {} };
    passwordRef = React.createRef();
    passwordConfirmRef = React.createRef();

    isValid = () => {
        var compulsoryFields = ['password', 'passwordConfirm'];
        var errors = {};
        const pw = this.passwordRef.current ? this.passwordRef.current.value : '';
        const pw2 = this.passwordConfirmRef.current ? this.passwordConfirmRef.current.value : '';
        if (pw !== pw2) {
            errors.passwordConfirm = "The passwords don't match";
        }
        compulsoryFields.forEach(function (field) {
            const node = field === 'password' ? this.passwordRef.current : this.passwordConfirmRef.current;
	    var value = (node && node.value ? node.value : '').trim();
            if (!value) {
                errors[field] = 'This field is required';
            }
        }.bind(this));
        this.setState({errors: errors});
		return Object.keys(errors).length === 0;
    };

    getFormData = () => {
        var data = {
		password: this.passwordRef.current ? this.passwordRef.current.value : '',
		passwordConfirm: this.passwordConfirmRef.current ? this.passwordConfirmRef.current.value : ''
        };
        return data;
    };

    render() {
        return (
			<div className="form-horizontal">
				{this.renderPassword('password', 'Password')}
				{this.renderPassword('passwordConfirm', 'Confirm Password')}
			</div>);
    }

    renderPassword(id, label) {
        return this.renderField(id, label,
            <input type="password" className="form-control" id={id} ref={id === 'password' ? this.passwordRef : this.passwordConfirmRef}/>
        );
    }

    renderField(id, label, field) {
        return (
			<div className={$c('form-group', {'has-error': id in this.state.errors})}>
				<label htmlFor={id} className="col-sm-4 control-label">{label}</label>
				<div className="col-sm-6">
					{field}
				</div>
			</div>);
    }
}

const ChangePassword = withRouter(ChangePasswordInner);
export { ChangePassword };
export default ChangePassword;
