/*global grecaptcha: false */
'use strict';

var React = require('react');
var content = require('./content');
var RawHTML = require('./RawHTML');
var countries = require('raw!../content/countries.txt').split("\n");
var $ = require('jquery');
var config  = require('./config');
var {Grid, Row, Col, Button} = require('react-bootstrap');
var {Navigation} = require('react-router');

var Role = {
    ROLE_DATA_PROVIDER: 12,
    options: [
        // [ id, dropdown text, community page text ]
        [1, "I am a patient",                      "Patient"],
        [2, "I am a proband",                      "Proband"],
        [3, "I am a clinical lab director",        "Clinical Lab Director"],
        [4, "I am a member of a diagnostic lab",   "Diagnostic Lab Staff"],
        [5, "I am a principal investigator",       "Principal Investigator"],
        [6, "I am a researcher",                   "Researcher"],
        [7, "I lead an advocacy group",            "Advocacy Group Leader"],
        [8, "I am a member of an advocacy group",  "Advocacy Group Member"],
        [9, "I am a genetic counselor",            "Genetic Counselor"],
        [10, "I am a clinical geneticist",         "Clinical Geneticist"],
        [11, "I am a clinician",                   "Clinician"],
        [12, "I represent a Data Provider",        "Data Provider"],
        [0, "Other"]
    ],
    other: id => parseInt(id) === 0 || parseInt(id) === 11,
    get: function(id) { return this.options.find(role => role[0] === parseInt(id)); }
};

var Signup = React.createClass({
    mixins: [Navigation],
    getInitialState: function () {
        return {
            submitted: null,
            success: null
        };
    },
    render: function () {
        var message;
        if (this.state.error != null) {
            message = (
				<div className="alert alert-danger">
					<p>{this.state.error}</p>
				</div>);
        }
        return (
            <Grid id="main-grid">
                <Row>
                    <Col sm={10} smOffset={1}  className="alert alert-warning">
                        <RawHTML ref='content' html={content.pages.signupMessage}/>
                    </Col>
                </Row>
                <Row id="message">
                    {message}
                </Row>
                <Row id="form">
                    <Col md={8} mdOffset={2}>
                        <SignupForm ref="contactForm"/>
                    </Col>
                </Row>
                <Row id="submit">
                    <Col md={6} mdOffset={3}>
                        <Button type="button" className="btn btn-primary btn-block" onClick={this.handleSubmit}>
                            Submit
                        </Button>
                    </Col>
                </Row>
            </Grid>);
    },

    handleChange: function (field, e) {
        var nextState = {};
        nextState[field] = e.target.checked;
        this.setState(nextState);
    },

    handleSubmit: function () {
        var self = this;
        var showSuccess = () => {this.transitionTo('/community', null, {registrationSuccess: true});};
        var showFailure = msg => {this.setState({error: msg});};

        var withGoogleMaps = function () {
            var geo = new google.maps.Geocoder();
            var formData = self.refs.contactForm.getFormData();
            var address = "" + formData.institution + "," + formData.city + "," + formData.state + "," + formData.country;

            var submit = function () {
                self.setState({submitted: formData});
                var url = config.backend_url + '/accounts/register/';

                var fd = new FormData();
                $.each(formData, function (k, v) {
                    fd.append(k, v);
                });

                var xhr = new XMLHttpRequest();
                xhr.onload = function () {
                    var responseData = JSON.parse(this.response);

                    if (this.status === 200 && responseData.success === true) {
                        showSuccess();
                    } else {
                        var message = responseData.error;
                        if (message === null) {
                            message = "Could not complete registration";
                        }
                        showFailure(message);
                    }
                };
                xhr.open('post', url);
                xhr.send(fd);
            };

            if (address.length > 3) {
                geo.geocode({address: address}, (results, status) => {
                    var loc;
                    if (status === google.maps.GeocoderStatus.OK) {
                        loc = results[0].geometry.location;
                        formData.latitude = loc.lat().toString();
                        formData.longitude = loc.lng().toString();
                    } else if (status === google.maps.GeocoderStatus.ZERO_RESULTS) {
                        showFailure("Please check your location information, or leave it blank.");
                        return;
                    } else if (status === google.maps.GeocoderStatus.OVER_QUERY_LIMIT) {
                        showFailure("Error checking your location information, please submit again.");
                        return;
                    }
                    submit();
                });
            } else {
                submit();
            }
        };

        if (this.refs.contactForm.isValid()) {
            google.load('maps', '3', {callback: withGoogleMaps, "other_params": "key=" + config.maps_key});
        } else {
            this.setState({error: "Some information was missing"});
        }
    }
});

function $c(staticClassName, conditionalClassNames) {
    var classNames = [];
    if (typeof conditionalClassNames === 'undefined') {
        conditionalClassNames = staticClassName;
    }
    else {
        classNames.push(staticClassName);
    }
    for (var className in conditionalClassNames) {
        if (conditionalClassNames[className]) {
            classNames.push(className);
        }
    }
    return classNames.join(' ');
}

var SignupForm = React.createClass({
    getInitialState: function () {
        return {errors: {}, file: '', imagePreviewUrl: null, captcha: "", otherRole: false};
    },
    componentDidMount: function () {
        var me = this;
        window.onRecaptchaLoad(function () {
            grecaptcha.render(me.refs.signupCAPTCHA.getDOMNode(), {sitekey: config.captcha_key, callback: function(resp) {
                me.setState({captcha: resp});
            }});
        });
    },
    isValid: function () {
        var errors = {};
        if (this.refs.role.getDOMNode().value === "NONE") {
            errors["role"] = "Please select a roll"; //eslint-disable-line dot-notation
        }
        if (this.refs.email.getDOMNode().value !== this.refs.email_confirm.getDOMNode().value) {
            errors["email_confirm"] = "The emails don't match"; //eslint-disable-line dot-notation
        }
        if (this.refs.password.getDOMNode().value !== this.refs.password_confirm.getDOMNode().value) {
            errors["password_confirm"] = "The passwords don't match"; //eslint-disable-line dot-notation
        }
        if (this.state.captcha === "") {
            errors["captcha"] = "No CAPTCHA entered"; //eslint-disable-line dot-notation
        }
        this.getCompulsoryFields().forEach(function (field) {
            var value = this.refs[field].getDOMNode().value.trim();
            if (!value) {
                errors[field] = 'This field is required';
            }
        }.bind(this));
        this.setState({errors: errors});

		return Object.keys(errors).length === 0;
    },
    getCompulsoryFields: function () {
        var fields = ['email', 'email_confirm', 'password', 'password_confirm', 'role'];
        if (!this.refs || !this.refs.role || parseInt(this.refs.role.getDOMNode().value) !== Role.ROLE_DATA_PROVIDER) {
            fields.push('firstName', 'lastName');
        }
        if (this.state.otherRole) {
            fields.push('role_other');
        }
        return fields;
    },
    getFormData: function () {
        var title = this.refs.titlemd.getDOMNode().checked && this.refs.titlemd.getDOMNode().value ||
            this.refs.titlephd.getDOMNode().checked && this.refs.titlephd.getDOMNode().value ||
            this.refs.titleother.getDOMNode().checked && this.refs.titlecustom.getDOMNode().value;

        var data = {
            "image": this.state.file,
            "email": this.refs.email.getDOMNode().value,
            "email_confirm": this.refs.email_confirm.getDOMNode().value,
            "password": this.refs.password.getDOMNode().value,
            "password_confirm": this.refs.password_confirm.getDOMNode().value,
            "firstName": this.refs.firstName.getDOMNode().value,
            "lastName": this.refs.lastName.getDOMNode().value,
            "title": title,
            "role": this.refs.role.getDOMNode().value,
            "role_other": this.state.otherRole ? this.refs.role_other.getDOMNode().value : Role.get(this.refs.role.getDOMNode().value)[2],
            "institution": this.refs.institution.getDOMNode().value,
            "city": this.refs.city.getDOMNode().value,
            "state": this.refs.state.getDOMNode().value,
            "country": this.refs.country.getDOMNode().value,
            "phone_number": this.refs.phone_number.getDOMNode().value,
            "hide_number": this.refs.hide_number.getDOMNode().checked,
            "hide_email": this.refs.hide_email.getDOMNode().checked,
            "captcha": this.state.captcha
        };
        return data;
    },
    handleImageChange(e) {
        e.preventDefault();

        let reader = new FileReader();
        let file = e.target.files[0];
        reader.onloadend = () => {
            if (file.size <= 4 * 1024 * 1024) {
                this.setState({
                    file: file,
                    imagePreviewUrl: reader.result,
                    imageTooBig: false
                });
            } else {
                this.setState({
                    file: null,
                    imagePreviewUrl: null,
                    imageTooBig: true
                });
            }
        };
        reader.readAsDataURL(file);
    },
    render: function () {
        var onChange = function() {
            var value = this.refs.role.getDOMNode().value;
            this.setState({otherRole: Role.other(value)});
        };
        return (
            <div className="form-horizontal" onChange={onChange.bind(this)}>
                {this.renderImageUpload('image', 'Profile picture')}
                {this.renderTextInput('email', 'Email')}
                {this.renderTextInput('email_confirm', 'Confirm Email')}
                {this.renderPassword('password', 'Password')}
                {this.renderPassword('password_confirm', 'Confirm Password')}
                {this.renderTextInput('firstName', 'First Name')}
                {this.renderTextInput('lastName', 'Last Name')}
                {this.renderRadioInlines('title', '', {
                    values: [{name: 'M.D.', ref: 'md'}, {name: 'Ph.D', ref: 'phd'}, {name: 'Other', ref: 'other'}]
                    , defaultCheckedValue: 'M.D.'
                })}
                {this.renderRoles()}
                {this.state.otherRole &&
                    <div className="slide-fade-in">{this.renderTextInput('role_other', <span style={{color: "#D00000"}}>Please Specify:</span>)}</div>}
                {this.renderTextInput('institution', 'Institution, Hospital or Company')}
                {this.renderTextInput('city', 'City')}
                {this.renderTextInput('state', 'State or Province')}
                {this.renderSelect('country', 'Country', countries.map(v => [v, v]))}
                {this.renderTextInput('phone_number', 'Phone number')}
                {this.renderCheckBox('hide_number', 'Hide my phone number on this website')}
                {this.renderCheckBox('hide_email', 'Hide my email address on this website')}
                {this.renderCAPTCHA('captcha', 'CAPTCHA *')}
			</div>);
    },
    renderImageUpload: function (id, label) {
        var {imagePreviewUrl, imageTooBig} = this.state;
        var imagePreview = null;
        var error = null;
        if (imagePreviewUrl) {
            imagePreview = (<img src={imagePreviewUrl} className="img-thumbnail" style={{maxHeight: '160px', maxWidth: '160px'}} />);
        }
        if (imageTooBig) {
            error = <p className="bg-danger">Please choose an image less than 4MB</p>;
        }
        return this.renderField(id, label,
            <div>
                <input onChange={this.handleImageChange} type="file" accept="image/*"/>
                {imagePreview}
                {error}
            </div>);
    },
    renderTextInput: function (id, label) {
        return this.renderField(id, label,
            <input type="text" className="form-control" id={id} ref={id}/>
        );
    },
    renderPassword: function (id, label) {
        return this.renderField(id, label,
            <input type="password" className="form-control" id={id} ref={id}/>
        );
    },
    renderTextarea: function (id, label) {
        return this.renderField(id, label,
            <textarea className="form-control" id={id} ref={id}/>
        );
    },
    renderSelect: function(id, label, opts) {
        var options = opts.map(value => <option key={id + value[0]} value={value[0]}>{value[1]}</option>);
        return this.renderField(id, label,
            <select className="form-control" id={id} ref={id}>
                <option key={id + "NONE"} value=""></option>
                {options}
            </select>
        );
    },
    renderRoles: function () {
        var id = 'role';
        var options = Role.options.map(value => <option key={id + value[0]} value={value[0]}>{value[1]}</option>);
        return this.renderField(id, 'Role',
            <select className="form-control" id={id} ref={id}>
                <option key={id + "NONE"} value="NONE">Choose one:</option>
                {options}
            </select>
        );
    },
    renderRadioInlines: function (id, label, kwargs) {
        var options = kwargs.values.map(function (value) {
            var defaultChecked = (value.name === kwargs.defaultCheckedValue);
            return (
				<label className="radio-inline">
					<input type="radio" ref={id + value.ref} name={id} value={value.name} defaultChecked={defaultChecked}/>
					{value.name}
				</label>);
        });
        options = <span className="col-sm-9">{options}</span>;
        var other =
            (<span className="col-sm-3">
            <input className="form-control" type="text" ref="titlecustom" name="titlecustom"/>
            </span>);
        var optionsWithOther = {options, other};
        return this.renderField(id, label, optionsWithOther);
    },
    renderCheckBox: function (id, label, defaultChecked = false) {
        var checkbox = (<label className="radio-inline">
            <input type='checkbox' ref={id} defaultChecked={defaultChecked}/>
            {label}
        </label>);
        return this.renderField(id, "", checkbox);
    },
    renderCAPTCHA: function(id, label) {
        return this.renderField(id, label, <div ref="signupCAPTCHA"></div>);
    },
    renderField: function (id, label, field) {
        return (
			<div className={$c('form-group', {'has-error': id in this.state.errors, 'required': this.getCompulsoryFields().includes(id)})}>
				<label htmlFor={id} className="col-sm-4 control-label">{label}</label>
				<div className="col-sm-6">
					{field}
            </div>
        </div>);
    }
});

module.exports = ({
    Signup: Signup,
    Role: Role,
    $c: $c
});
