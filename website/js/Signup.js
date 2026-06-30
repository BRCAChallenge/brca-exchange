'use strict';

import React from 'react';
import content from './content';
import RawHTML from './RawHTML';
var countries = require('raw-loader!../content/countries.txt')["default"].split("\n");
var $ = require('jquery');
var _ = require('underscore');
import config from './config';
const { Container: Grid, Row, Col, Button} = require('react-bootstrap');
const { withRouter } = require('react-router-dom');

// ---- Google Maps loader (async, no jsapi) ----
function loadGoogleMaps(key, cb, onNoKey) {
  if (!key) {
    // If no key, keep registration functional (just skip lat/lng enrichment)
    if (typeof onNoKey === 'function') onNoKey();
    return;
  }
  if (window.google && window.google.maps) {
    cb();
    return;
  }
  const script = document.createElement('script');
  window.__signupInitMap = cb;
  script.src = `https://maps.googleapis.com/maps/api/js?key=${encodeURIComponent(
    key
  )}&callback=__signupInitMap&loading=async`;
  script.async = true;
  script.defer = true;
  document.head.appendChild(script);
}

export const Role = {
    ROLE_DATA_PROVIDER: 12,
    options: [
        // [ id, dropdown text, community page text ]
        [1, "I am a concerned member of the public",    "Concerned Person"],
        [3, "I am a clinical lab director",             "Clinical Lab Director"],
        [4, "I am a member of a diagnostic lab",        "Diagnostic Lab Staff"],
        [5, "I am a principal investigator",            "Principal Investigator"],
        [6, "I am a researcher",                        "Researcher"],
        [7, "I lead an advocacy group",                 "Advocacy Group Leader"],
        [8, "I am a member of an advocacy group",       "Advocacy Group Member"],
        [9, "I am a genetic counselor",                 "Genetic Counselor"],
        [10, "I am a clinical geneticist",              "Clinical Geneticist"],
        [11, "I am a clinician",                        "Clinician"],
        [12, "I represent a Data Provider",             "Data Provider"],
        [0, "Other"]
    ],
    other: id => parseInt(id) === 0,
    get: function(id) { return this.options.find(role => role[0] === parseInt(id)); }
};

class SignupInner extends React.Component {
    contentRef = React.createRef();
    contactFormRef = React.createRef();
    constructor(props) {
        super(props);
        this.contactFormRef = React.createRef();
    }

    state = {
        submitted: null,
        success: null
    };

    render() {
        var message;
        if (this.state.error != null) {
            let fieldErrors = _.map(this.state.fieldErrors, (err, i) => (<li key={`fe-${i}`}>{err}</li>));
            fieldErrors = _.values(_.groupBy(fieldErrors, (item, index) => Math.floor(index / 2) )).map((group, i) => <Col key={`fec-${i}`} md={3}><ul>{group}</ul></Col>);
            message = (
				<div className="alert alert-danger">
					<Row><Col md={6}>{this.state.error}</Col></Row>
                    <Row><Col md={1} />{fieldErrors}</Row>
				</div>);
            window.scrollTo(0, 0);
        }
        return (
            <Grid id="main-grid">
                <Row>
                    <Col sm={{ span: 10, offset: 1 }}  className="alert alert-warning">
                        <RawHTML ref={this.contentRef} html={content.pages.signupMessage}/>
                    </Col>
                </Row>
                <Row id="message">
                    {message}
                </Row>
                <Row id="form">
                    <Col md={{ span: 8, offset: 2 }}>
                        <SignupForm ref={this.contactFormRef}/>
                    </Col>
                </Row>
                <Row id="submit">
                    <Col md={{ span: 6, offset: 3 }}>
                        <Button type="button" className="btn btn-primary btn-block" onClick={this.handleSubmit}>
                            Submit
                        </Button>
                    </Col>
                </Row>
            </Grid>);
    }

    handleChange = (field, e) => {
        var nextState = {};
        nextState[field] = e.target.checked;
        this.setState(nextState);
    };

    handleSubmit = () => {
        var showSuccess = () => {
		this.props.history.push({ pathname: '/community', search: '?registrationSuccess=true' });
	};
        var showFailure = msg => {this.setState({error: msg});};
        var address = '';


        var submit = function (formData) {
            this.setState({submitted: formData});
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
        }.bind(this);

        var addLeadingCommaIfNecessary = function(address) {
            if (address.length > 0) {
                address += ", ";
            }
            return address;
        };

        var determineAddressFromCityStateCountry = function(formData) {
            if (formData.city) {
                address = addLeadingCommaIfNecessary(address);
                address += formData.city;
            }
            if (formData.state) {
                address = addLeadingCommaIfNecessary(address);
                address += formData.state;
            }
            if (formData.country) {
                address = addLeadingCommaIfNecessary(address);
                address += formData.country;
            }
            return address;
        };

        var getLatLng = function(address, formData) {
            var geo = new window.google.maps.Geocoder();
            geo.geocode({address: address}, (results, status) => {
                var loc;
                if (status === window.google.maps.GeocoderStatus.OK) {
                    loc = results[0].geometry.location;
                    formData.latitude = loc.lat().toString();
                    formData.longitude = loc.lng().toString();
                } else {
                    console.log("Error parsing address.");
                }
                submit(formData);
            });
        };

        var withGoogleMaps = function () {
            var formData = this.contactFormRef.current.getFormData();

            address = determineAddressFromCityStateCountry(formData);
            if (address.length > 3) {
                getLatLng(address, formData);
            } else if (formData.institution.length > 3) {
                getLatLng(formData.institution, formData);
            } else {
                submit(formData);
            }
        }.bind(this);

        var formErrors = this.contactFormRef.current.getFormErrors();
        if (formErrors === false) {
	    loadGoogleMaps(config.maps_key, withGoogleMaps, () => submit(this.contactFormRef.current.getFormData()));
        } else {
            this.setState({error: <strong>Some information was missing:</strong>, fieldErrors: formErrors });
        }
    };
}

export function $c(staticClassName, conditionalClassNames) {
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

class SignupForm extends React.Component {
    state = {errors: {}, file: '', imagePreviewUrl: null, captcha: "", otherRole: false};
    _refs = {};
    setRef = (name) => (el) => {
        if (el) this._refs[name] = el;
        else delete this._refs[name];
    };

    componentDidMount() {
        const renderCaptcha = () => {
            if (!this._refs.signupCAPTCHA || !window.grecaptcha) return;
            window.grecaptcha.render(this._refs.signupCAPTCHA, {
                sitekey: config.captcha_key,
                callback: (resp) => this.setState({ captcha: resp })
            });
        };
        if (typeof window.onRecaptchaLoad === 'function') {
            window.onRecaptchaLoad(renderCaptcha);
        } else {
            // if your recaptcha script loads before this component mounts
            renderCaptcha();
        }
    }

    getFormErrors() {
        var errors = {};
        if (!this._refs.role || this._refs.role.value === "NONE") {
            errors.role = <span>Please select a <strong>Role</strong></span>;
        }
        if ((this._refs.email?.value || "") !== (this._refs.email_confirm?.value || "")) {
            errors["email_confirm"] = <span>The <strong>emails</strong> don't match</span>;
        }
        if ((this._refs.password?.value || "") !== (this._refs.password_confirm?.value || "")) {
            errors["password_confirm"] = <span>The <strong>passwords</strong> don't match</span>;
        }
        if (this.state.captcha === "") {
            errors.captcha = <span>No <strong>CAPTCHA</strong> entered</span>;
        }
        this.getCompulsoryFields().forEach(function (field) {
            var value = (this._refs[field]?.value || "").trim();
            if (!value) {
                errors[field] = <span><strong>{ field.replace(/([A-Z])/g, ' $1').replace(/^./, function(str) { return str.toUpperCase(); }) }</strong> is required</span>;
            }
        }.bind(this));
        this.setState({errors: errors});

		if (Object.keys(errors).length === 0) {
            return false;
        } else {
            return errors;
        }
    }

    getCompulsoryFields() {
        var fields = ['email', 'password', 'role'];
	const roleVal = this._refs.role ? parseInt(this._refs.role.value) : NaN;
        if (!this._refs.role || roleVal !== Role.ROLE_DATA_PROVIDER) {
            fields.push('firstName', 'lastName');
        }
        if (this.state.otherRole) {
            fields.push('role_other');
        }
        return fields;
    }

    getFormData() {
        var title =
            (this._refs.titlemd?.checked && this._refs.titlemd?.value) ||
            (this._refs.titlephd?.checked && this._refs.titlephd?.value) ||
            (this._refs.titleother?.checked && this._refs.titlecustom?.value);

        const roleVal = this._refs.role ? this._refs.role.value : "NONE";
        const roleLabel = Role.get(roleVal) ? Role.get(roleVal)[2] : "";

	var data = {
            "image": this.state.file,
            "email": this._refs.email?.value,
            "email_confirm": this._refs.email_confirm?.value,
            "password": this._refs.password?.value,
            "password_confirm": this._refs.password_confirm?.value,
            "firstName": this._refs.firstName?.value,
            "lastName": this._refs.lastName?.value,
	    "title": title,
            "role": roleVal,
            "role_other": this.state.otherRole ? (this._refs.role_other?.value) : roleLabel,
            "institution": this._refs.institution?.value,
            "city": this._refs.city?.value,
            "state": this._refs.state?.value,
            "country": this._refs.country?.value,
            "phone_number": this._refs.phone_number?.value,
            "hide_number": !!this._refs.hide_number?.checked,
            "hide_email": !!this._refs.hide_email?.checked,
	    "captcha": this.state.captcha
        };
        return data;
    }

    handleImageChange = (e) => {
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
    };

    render() {
        var onChange = function() {
            var value = this._refs.role ? this._refs.role.value : "NONE";
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
    }
    renderImageUpload(id, label) {
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
    }
    renderTextInput(id, label) {
        return this.renderField(id, label,
            <input type="text" className="form-control" id={id} ref={this.setRef(id)}/>
        );
    }
    renderPassword(id, label) {
        return this.renderField(id, label,
            <input type="password" className="form-control" id={id} ref={this.setRef(id)}/>
        );
    }
    renderTextarea(id, label) {
        return this.renderField(id, label,
            <textarea className="form-control" id={id} ref={this.setRef(id)}/>
        );
    }
    renderSelect(id, label, opts) {
        var options = opts.map(value => <option key={id + value[0]} value={value[0]}>{value[1]}</option>);
        return this.renderField(id, label,
            <select className="form-control" id={id} ref={this.setRef(id)}>
                <option key={id + "NONE"} value="" />
                {options}
            </select>
        );
    }
    renderRoles() {
        var id = 'role';
        var options = Role.options.map(value => <option key={id + value[0]} value={value[0]}>{value[1]}</option>);
        return this.renderField(id, 'Role',
            <select className="form-control" id={id} ref={this.setRef(id)}>
                <option key={id + "NONE"} value="NONE">Choose one:</option>
                {options}
            </select>
        );
    }
    renderRadioInlines(id, label, kwargs) {
        var options = kwargs.values.map((value) => {
            var defaultChecked = (value.name === kwargs.defaultCheckedValue);
            if (value.name === "Other") {
                return (
                    <label key={`${id}-${value.ref}`} className="radio-inline">
                        <input type="radio" ref={this.setRef(id + value.ref)} name={id} value={value.name} defaultChecked={defaultChecked}/>
                        {value.name + ':'}<input className="form-control" type="text" ref={this.setRef("titlecustom")} name="titlecustom"/>
                    </label>
                );
            } else {
                return (<label key={`${id}-${value.ref}`} className="radio-inline">
                    <input type="radio" ref={this.setRef(id + value.ref)} name={id} value={value.name} defaultChecked={defaultChecked}/>
                        {value.name}
                    </label>
                );
            }
        });
        options = (<span className="col-xs-12">{options}</span>);
        return this.renderField(id, label, options);
    }
    renderCheckBox(id, label, defaultChecked = false) {
        var checkbox = (<label className="radio-inline">
            <input type='checkbox' ref={this.setRef(id)} defaultChecked={defaultChecked}/>
            {label}
        </label>);
        return this.renderField(id, "", checkbox);
    }
    renderCAPTCHA(id, label) {
        return this.renderField(id, label, <div ref={this.setRef("signupCAPTCHA")} />);
    }
    renderField(id, label, field) {
        return (
			<div className={$c('form-group', {'has-error': id in this.state.errors, 'required': this.getCompulsoryFields().includes(id)})}>
				<label htmlFor={id} className="col-sm-4 control-label">{label}</label>
				<div className="col-sm-6">
					{field}
                </div>
            </div>
        );
    }
}

const Signup = withRouter(SignupInner);
export { Signup };
export default Signup;
