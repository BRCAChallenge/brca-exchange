'use strict';

import React from 'react';
import ReactDOM from 'react-dom';
import $ from 'jquery';
import _ from 'underscore';
import config from './config';
import { Container as Grid, Row, Col, Button } from 'react-bootstrap';
import { withRouter } from 'react-router-dom';
import auth from './auth';
import { Role, $c } from './Signup';
var countries = require('raw-loader!../content/countries.txt').default.split('\n');

class ProfileInner extends React.Component {
    state = { success: null, error: null, fieldErrors: null };
    contactFormRef = React.createRef();

    componentDidMount() {
        if (!auth.loggedIn()) {
            const target = (this.props.location && this.props.location.pathname) || '/profile';
            this.props.history.replace({
                pathname: '/signin',
                search: `?target=${encodeURIComponent(target)}`
            });
        }
    }

    render() {
        var message;
        if (this.state.error != null) {
            let fieldErrors = _.map(this.state.fieldErrors, err => (<li>{err}</li>));
            fieldErrors = _.values(_.groupBy(fieldErrors, (item, index) => Math.floor(index / 2) )).map(group => <Col md={3}><ul>{group}</ul></Col>);
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
                    <div className='text-center Variant-detail-title'>
                        <h3>Update your profile</h3>
                    </div>
                </Row>
                <Row id="message">
                    {message}
                </Row>
                <Row id="form">
                    <Col md={{ span: 8, offset: 2 }}>
                        <EditProfileForm ref={this.contactFormRef} history={this.props.history}/>
                    </Col>
                </Row>
                <Row id="submit">
                    <Col md={{ span: 6, offset: 3 }}>
                        <Button type="button" className="btn btn-primary btn-block" onClick={this.handleSubmit}>
                            Update
                        </Button>
                    </Col>
                </Row>
            </Grid>);
    }

    handleChange = (field, e)  => {
        var nextState = {};
        nextState[field] = e.target.checked;
        this.setState(nextState);
    };

    handleSubmit = () => {
        var self = this;
	const form = this.contactFormRef.current;
        const subscribe = form ? form.getSubscribeAction() : undefined;
        var showSuccess = () => {
            const qs = new URLSearchParams();
            qs.set('updateSuccess', 'true');
            if (subscribe !== undefined) qs.set('subscribe', String(!!subscribe));
            this.props.history.push({ pathname: '/community', search: `?${qs.toString()}` });
        };
        var showFailure = msg => {this.setState({error: msg || "An error occurred."});};

        var withGoogleMaps = function () {
            var geo = new google.maps.Geocoder();
            var formData = form.getFormData();
            var address = "" + formData.institution + "," + formData.city + "," + formData.state + "," + formData.country;

            var submit = function() {
                self.setState({submitted: formData});
                var url = config.backend_url + '/accounts/update/';

                var fd = new FormData();
                $.each(formData, function (k, v) {
                    fd.append(k, v);
                });

                if(subscribe !== undefined) {
                    fd.append("subscribe", subscribe);
                }

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
                xhr.setRequestHeader('Authorization', 'Bearer ' + auth.token());
                xhr.send(fd);
            };

            if (address.length > 3) {
                geo.geocode({address: address}, (results, status) => {
                    var loc;
                    if (status === google.maps.GeocoderStatus.OK) {
                        loc = results[0].geometry.location;
                        formData.latitude = loc.lat().toString();
                        formData.longitude = loc.lng().toString();
                    }
                    /* else if (status === google.maps.GeocoderStatus.ZERO_RESULTS) {
                        showFailure("Please check your location information, or leave it blank.");
                        return;
                    } else if (status === google.maps.GeocoderStatus.OVER_QUERY_LIMIT) {
                        showFailure("Error checking your location information, please submit again.");
                        return;
                    } */
                    submit();
                });
            } else {
                submit();
            }
        };

        var formErrors = form && form.getFormErrors ? form.getFormErrors() : { form: 'Form not ready' };
        if (formErrors === false) {
            google.load('maps', '3', {callback: withGoogleMaps, "other_params": "key=" + config.maps_key});
        } else {
            this.setState({error: <strong>Some information was missing:</strong>, fieldErrors: formErrors });
        }
    };
}

class EditProfileForm extends React.Component {
    state = { errors: {}, data: {}, imagePreviewUrl: '', otherRole: false };
    _refs = {};
    setRef = (name) => (el) => { if (el) this._refs[name] = el; else delete this._refs[name]; };
    getNode = (name) => this._refs[name] || null;

    componentDidMount() {
        this.retrieveProfile();
    }

    getSubscribeAction() {
        // returns undefined to do nothing, true to subscribe, false to unsubscribe
        return (
            this.state.mailingList !== this.state.oldMailingList
                ? this.state.mailingList
                : undefined
        );
    }

    retrieveProfile() {
        var url = config.backend_url + '/accounts/get/';
        var token = auth.token();
        var tokenValue = 'Bearer ' + token;
        var saveProfileData = (data) => {
            var imagePreviewUrl = '';
            if (data.user.has_image) {
                imagePreviewUrl = config.backend_url + '/site_media/media/' + data.user.id;
            }
            var otherRole = Role.other(data.user.role);
            this.setState({data: data.user, mailingList: data.mailinglist, oldMailingList: data.mailinglist, imagePreviewUrl: imagePreviewUrl, otherRole: otherRole});
        };
        $.ajax({
            type: 'GET',
            headers: {'Authorization': tokenValue},
            url: url,
            success: function (data) {
                saveProfileData(data);
            },
            error: function () {
                if (this.props.history) {
                    this.props.history.replace({ pathname: '/signin', search: `?target=${encodeURIComponent('/profile')}` });
                }
            }.bind(this)
        });
    }

    getFormErrors() {
        var errors = {};
	const roleNode = this.getNode('role');
        if (roleNode && roleNode.value === "NONE") {
            errors.role = <span>Please select a <strong>Role</strong></span>;
        }
        if (this.state.captcha === "") {
            errors.captcha = <span>No <strong>CAPTCHA</strong> entered</span>;
        }
        this.getCompulsoryFields().forEach(function (field) {
            const node = this.getNode(field);
            var value = (node && node.value ? node.value : '').trim();
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
        var fields = [];
	const roleNode = this.getNode('role');
        if (!roleNode || parseInt(roleNode.value) !== Role.ROLE_DATA_PROVIDER) {
            fields.push('firstName', 'lastName');
        }
        if (this.state.otherRole) {
            fields.push('role_other');
        }
        return fields;
    }

    getFormData() {
        const title = 
	    (this.getNode('titlemd') && this.getNode('titlemd').checked && this.getNode('titlemd').value) ||
            (this.getNode('titlephd') && this.getNode('titlephd').checked && this.getNode('titlephd').value) ||
            (this.getNode('titleother') && this.getNode('titleother').checked && (this.getNode('titlecustom') ? this.getNode('titlecustom').value : ''));

        var data = {
            "image": this.state.file,
            "deleteImage": this.state.imageDelete,
            "firstName": (this.getNode('firstName') && this.getNode('firstName').value) || '',
            "lastName": (this.getNode('lastName') && this.getNode('lastName').value) || '',
            "title": title,
            "role": (this.getNode('role') && this.getNode('role').value) || '',
	    "role_other": this.state.otherRole
                ? ((this.getNode('role_other') && this.getNode('role_other').value) || '')
                : (Role.get(((this.getNode('role') && this.getNode('role').value) || 0)) || [])[2],
            "institution": (this.getNode('institution') && this.getNode('institution').value) || '',
            "city": (this.getNode('city') && this.getNode('city').value) || '',
            "state": (this.getNode('state') && this.getNode('state').value) || '',
            "country": (this.getNode('country') && this.getNode('country').value) || '',
            "phone_number": (this.getNode('phone_number') && this.getNode('phone_number').value) || '',
            "hide_number": !!(this.getNode('hide_number') && this.getNode('hide_number').checked),
            "hide_email": !!(this.getNode('hide_email') && this.getNode('hide_email').checked)
        };
        return data;
    }

    handleImageChange(e) {
        e.preventDefault();

        let reader = new FileReader();
        let file = e.target.files[0];
        reader.onloadend = () => {
            if (file.size <= 4 * 1024 * 1024) {
                this.setState({
                    file: file,
                    imagePreviewUrl: reader.result,
                    imageTooBig: false,
                    imageDelete: null
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
    }
    render() {
        return (
        <div className="form-horizontal">
            {this.renderImageUpload('image', 'Profile picture')}
            {
                /* {this.renderPassword('password', 'Password')}
                   {this.renderPassword('password_confirm', 'Confirm Password')}
                */
            }
            {this.renderTextInput('firstName', 'First Name', this.state.data.firstName)}
            {this.renderTextInput('lastName', 'Last Name', this.state.data.lastName)}
            {this.renderRadioInlines('title', '', {
                values: [{name: 'M.D.', ref: 'md'}, {name: 'Ph.D', ref: 'phd'}, {name: 'Other', ref: 'other'}]
                , defaultCheckedValue: this.state.data.title
            })}
            {this.renderRoles(this.state.data.role)}
            {this.state.otherRole &&
                <div className="slide-fade-in">{this.renderTextInput('role_other', <span style={{color: "#D00000"}}>Please Specify:</span>, this.state.data.role_other)}</div>}
            {this.renderTextInput('institution', 'Institution, Hospital or Company', this.state.data.institution)}
            {this.renderTextInput('city', 'City', this.state.data.city)}
            {this.renderTextInput('state', 'State or Province', this.state.data.state)}
            {this.renderSelect('country', 'Country', countries.map(v => [v, v]), this.state.data.country)}
            {this.renderTextInput('phone_number', 'Phone number', this.state.data.phone_number)}
            {this.renderCheckBox('hide_number', "Don't display my phone number on this website", this.state.data.hide_number)}
            {this.renderCheckBox('hide_email', "Don't display my email on this website", this.state.data.hide_email)}
            {this.renderMailingList('mailinglist', "Subscribed to mailing list?", this.state.mailingList)}
        </div>);
    }

    renderImageUpload(id, label) {
        var handleImageDelete = ()=>
            this.setState({
                imageDelete: true,
                imagePreviewUrl: '',
                file: null
            });
        var {imagePreviewUrl, imageTooBig} = this.state;
        var imagePreview = null;
        var error = null;
        if (imagePreviewUrl) {
            imagePreview = (
                <div>
                    <div><img src={imagePreviewUrl} className="img-thumbnail"
                              style={{maxHeight: '160px', maxWidth: '160px'}}/></div>
                    <div ><Button variant="link" onClick={handleImageDelete}>Remove picture</Button></div>
                </div>
            );
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

    renderTextInput(id, label, defaultValue) {
        var handleChange = (e) => {var oldData = this.state.data; oldData[id] = e.target.value; this.setState({data: oldData});};
        return this.renderField(id, label,
            <input type="text" className="form-control" id={id} ref={this.setRef(id)} value={defaultValue || ''} onChange={handleChange}/>
        );
    }
/*
    renderPassword: function (id, label) {
        return this.renderField(id, label,
            <input type="password" className="form-control" id={id} ref={id}/>
        );
    },
*/
    renderTextarea(id, label, defaultValue) {
        var handleChange = (e) => {var oldData = this.state.data; oldData[id] = e.target.value; this.setState({data: oldData});};
        return this.renderField(id, label,
            <textarea className="form-control" id={id} ref={this.setRef(id)} value={defaultValue} onChange={handleChange}/>
        );
    }
    renderRoles(defaultValue) {
        var id = 'role';
        var handleChange = (e) => {
            var oldData = this.state.data;
            oldData[id] = e.target.value;
            this.setState({otherRole: Role.other(e.target.value)});
        };
        var options = Role.options.map(value => <option key={id + value[0]} value={value[0]}>{value[1]}</option>);
        return this.renderField(id, 'Role',
            <select className="form-control" id={id} ref={this.setRef(id)} value={defaultValue || ''} onChange={handleChange}>
                {options}
            </select>
        );
    }

    renderSelect(id, label, opts, defaultValue) {
        var handleChange = (e) => {var oldData = this.state.data; oldData[id] = e.target.value; this.setState({data: oldData});};
        var options = opts.map(value => <option key={id + value[0]} value={value[0]}>{value[1]}</option>);
        return this.renderField(id, label,
            <select className="form-control" id={id} ref={this.setRef(id)} value={defaultValue || ''} onChange={handleChange}>
                <option key={id + "NONE"} value=""></option>
                {options}
            </select>
        );
    }

    renderRadioInlines(id, label, kwargs) {
        var handleTextChange = (e) => {var oldData = this.state.data; oldData[id] = e.target.value; this.setState({data: oldData});};
        var otherValue = kwargs.defaultCheckedValue;
		// XXX Not sure why eslint flags this bind, because 'this' is used in the handlers within the
		// body of the function.
        var options = kwargs.values.map((value) => { //eslint-disable-line no-extra-bind
            var handleRadioChange = () => {var oldData = this.state.data; oldData[id] = value.name; this.setState({data: oldData}); };
            var defaultChecked = false;
            if (value.name === kwargs.defaultCheckedValue) {
                defaultChecked = true;
                otherValue = '';
            }
            if (value.name === 'Other' && !kwargs.values.some(opt => opt.name === kwargs.defaultCheckedValue)) {defaultChecked = true;}
            return (
				<label className="radio-inline">
					<input type="radio" ref={this.setRef(id + value.ref)} name={id} value={value.name} checked={defaultChecked} onChange={handleRadioChange}/>
					{value.name}
				</label>);
        });
        options = <span className="col-sm-9">{options}</span>;
        var other =
            (<span className="col-sm-3">
            <input className="form-control" type="text" ref={this.setRef('titlecustom')} name="titlecustom" value={otherValue || ''} onChange={handleTextChange}/>
            </span>);
        return this.renderField(id, label, <>{options}{other}</>);
    }

    renderCheckBox(id, label, defaultValue) {
        var handleChange = (e) => {var oldData = this.state.data; oldData[id] = e.target.checked; this.setState({data: oldData});};
        var checkbox = (<label className="radio-inline">
            <input type='checkbox' ref={this.setRef(id)} checked={!!defaultValue} onChange={handleChange} />
            &nbsp;{label}
        </label>);
        return this.renderField(id, "", checkbox);
    }

    renderMailingList(id, label, defaultValue) {
        var handleChange = (e) => {this.setState({mailingList: !!e.target.checked});};
        var checkbox = (<label className="radio-inline">
            <input type='checkbox' ref={this.setRef(id)} checked={!!defaultValue} onChange={handleChange} />
            &nbsp;{label}
        </label>);
        return this.renderField(id, "", checkbox);
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

const Profile = withRouter(ProfileInner);
export { Profile };
export default Profile;
