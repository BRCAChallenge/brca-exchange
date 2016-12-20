'use strict';

var React = require('react');
var $ = require('jquery');
var _ = require('underscore');
var config  = require('./config');
var {Grid, Row, Col, Button} = require('react-bootstrap');
var {Navigation} = require('react-router');
var auth = require('./auth');
var {Role, $c} = require('./Signup');
var countries = require('raw!../content/countries.txt').split('\n');

var Profile = React.createClass({
    statics: {
        willTransitionTo: function (transition) {
            if (!auth.loggedIn()) {
                transition.redirect('/signin', {}, {
                    target: transition.path
                });
            }
        }
    },
    mixins: [Navigation],
    getInitialState: function () {
        return {
            success: null,
        };
    },
    render: function () {
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
                    <Col md={8} mdOffset={2}>
                        <EditProfileForm ref="contactForm"/>
                    </Col>
                </Row>
                <Row id="submit">
                    <Col md={6} mdOffset={3}>
                        <Button type="button" className="btn btn-primary btn-block" onClick={this.handleSubmit}>
                            Update
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
        var showSuccess = () => this.transitionTo('/community', null, {updateSuccess: true, subscribe: this.refs.contactForm.getSubscribeAction()});
        var showFailure = msg => {this.setState({error: msg || "An error occurred."});};

        var withGoogleMaps = function () {
            var geo = new google.maps.Geocoder();
            var formData = self.refs.contactForm.getFormData();
            var address = "" + formData.institution + "," + formData.city + "," + formData.state + "," + formData.country;

            var submit = function() {
                self.setState({submitted: formData});
                var url = config.backend_url + '/accounts/update/';

                var fd = new FormData();
                $.each(formData, function (k, v) {
                    fd.append(k, v);
                });

                if(self.refs.contactForm.getSubscribeAction() !== undefined) {
                    fd.append("subscribe", self.refs.contactForm.getSubscribeAction());
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
                xhr.setRequestHeader('Authorization', 'JWT ' + auth.token());
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

        var formErrors = this.refs.contactForm.getFormErrors();
        if (formErrors === false) {
            google.load('maps', '3', {callback: withGoogleMaps, "other_params": "key=" + config.maps_key});
        } else {
            this.setState({error: <strong>Some information was missing:</strong>, fieldErrors: formErrors });
        }
    }
});

var EditProfileForm = React.createClass({
    mixins: [Navigation],
    componentDidMount: function () {
        this.retrieveProfile();
    },
    getSubscribeAction: function () {
        // returns undefined to do nothing, true to subscribe, false to unsubscribe
        return (
            this.state.mailingList !== this.state.oldMailingList
                ? this.state.mailingList
                : undefined
        );
    },
    retrieveProfile: function () {
        var url = config.backend_url + '/accounts/get/';
        var token = auth.token();
        var tokenValue = 'JWT ' + token;
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
                this.transitionTo('/signin', {}, {target: '/profile'});
            }.bind(this)
        });
    },
    getInitialState: function () {
        return {errors: {}, data: {}};
    },
    getFormErrors: function () {
        var errors = {};
        if (this.refs.role.getDOMNode().value === "NONE") {
            errors.role = <span>Please select a <strong>Roll</strong></span>;
        }
        if (this.state.captcha === "") {
            errors.captcha = <span>No <strong>CAPTCHA</strong> entered</span>;
        }
        this.getCompulsoryFields().forEach(function (field) {
            var value = this.refs[field].getDOMNode().value.trim();
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
    },
    getCompulsoryFields: function () {
        var fields = [];
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
            "deleteImage": this.state.imageDelete,
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
            "hide_email": this.refs.hide_email.getDOMNode().checked
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
    },
    render: function () {
        var onChange = function() {
            var value = this.refs.role.getDOMNode().value;
            this.setState({otherRole: Role.other(value)});
        };
        return (
        <div className="form-horizontal" onChange={onChange.bind(this)}>
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
    },
    renderImageUpload: function (id, label) {
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
                    <div ><Button bsStyle="link" onClick={handleImageDelete}>Remove picture</Button></div>
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
    },
    renderTextInput: function (id, label, defaultValue) {
        var handleChange = () => {var oldData = this.state.data; oldData[id] = this.refs[id].getDOMNode().value; this.setState({data: oldData});};
        return this.renderField(id, label,
            <input type="text" className="form-control" id={id} ref={id} value={defaultValue} onChange={handleChange}/>
        );
    },
/*
    renderPassword: function (id, label) {
        return this.renderField(id, label,
            <input type="password" className="form-control" id={id} ref={id}/>
        );
    },
*/
    renderTextarea: function (id, label, defaultValue) {
        var handleChange = () => {var oldData = this.state.data; oldData[id] = this.refs[id].value; this.setState({data: oldData});};
        return this.renderField(id, label,
            <textarea className="form-control" id={id} ref={id} value={defaultValue} onChange={handleChange}/>
        );
    },
    renderRoles: function (defaultValue) {
        var id = 'role';
        var handleChange = () => {
            var oldData = this.state.data;
            oldData[id] = this.refs[id].value;
            this.setState({data: oldData});
        };
        var options = Role.options.map(value => <option key={id + value[0]} value={value[0]}>{value[1]}</option>);
        return this.renderField(id, 'Role',
            <select className="form-control" id={id} ref={id} value={defaultValue} onChange={handleChange}>
                {options}
            </select>
        );
    },
    renderSelect: function(id, label, opts, defaultValue) {
        var handleChange = () => {var oldData = this.state.data; oldData[id] = this.refs[id].checked; this.setState({data: oldData});};
        var options = opts.map(value => <option key={id + value[0]} value={value[0]}>{value[1]}</option>);
        return this.renderField(id, label,
            <select className="form-control" id={id} ref={id} value={defaultValue} onChange={handleChange}>
                <option key={id + "NONE"} value=""></option>
                {options}
            </select>
        );
    },
    renderRadioInlines: function (id, label, kwargs) {
        var handleTextChange = () => {var oldData = this.state.data; oldData[id] = this.refs.titlecustom.getDOMNode().value; this.setState({data: oldData});};
        var otherValue = kwargs.defaultCheckedValue;
		// XXX Not sure why eslint flags this bind, because 'this' is used in the handlers within the
		// body of the function.
        var options = kwargs.values.map(function (value) { //eslint-disable-line no-extra-bind
            var handleRadioChange = () => {var oldData = this.state.data; oldData[id] = value.name; this.setState({data: oldData}); };
            var defaultChecked = false;
            if (value.name === kwargs.defaultCheckedValue) {
                defaultChecked = true;
                otherValue = '';
            }
            if (value.name === 'Other' && !kwargs.values.some(opt => opt.name === kwargs.defaultCheckedValue)) {defaultChecked = true;}
            return (
				<label className="radio-inline">
					<input type="radio" ref={id + value.ref} name={id} value={value.name} checked={defaultChecked} onChange={handleRadioChange}/>
					{value.name}
				</label>);
        }.bind(this));
        options = <span className="col-sm-9">{options}</span>;
        var other =
            (<span className="col-sm-3">
            <input className="form-control" type="text" ref="titlecustom" name="titlecustom" value={otherValue} onChange={handleTextChange}/>
            </span>);
        var optionsWithOther = {options, other};
        return this.renderField(id, label, optionsWithOther);
    },
    renderCheckBox: function (id, label, defaultValue) {
        var handleChange = () => {var oldData = this.state.data; oldData[id] = this.refs[id].checked; this.setState({data: oldData});};
        var checkbox = (<label className="radio-inline">
            <input type='checkbox' ref={id} checked={defaultValue} onChange={handleChange} />
            &nbsp;{label}
        </label>);
        return this.renderField(id, "", checkbox);
    },
    renderMailingList: function (id, label, defaultValue) {
        var handleChange = () => {this.setState({mailingList: !!this.refs[id].getDOMNode().checked});};
        var checkbox = (<label className="radio-inline">
            <input type='checkbox' ref={id} checked={defaultValue} onChange={handleChange} />
            &nbsp;{label}
        </label>);
        return this.renderField(id, "", checkbox);
    },
    renderField: function (id, label, field) {
        return (
			<div className={$c('form-group', {'has-error': id in this.state.errors})}>
				<label htmlFor={id} className="col-sm-4 control-label">{label}</label>
				<div className="col-sm-6">
					{field}
				</div>
			</div>);
    }
});

module.exports = ({
    Profile: Profile
});
