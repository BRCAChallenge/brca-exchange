'use strict';

var React = require('react');
var backend = require('./backend');
var content = require('./content');
var RawHTML = require('./RawHTML');
var $ = require('jquery');
var config  = require('./config')
var {Grid, Row, Col, Button} = require('react-bootstrap');
var {Navigation} = require('react-router');
var auth = require('./auth');
var {Signup, AFFILIATION, trim, $c} = require('./Signup');

var Profile = React.createClass({
    statics: {
        willTransitionTo: function (transition, params, query) {
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
        }
    },
    render: function () {
        var message;
        if (this.state.error != null) {
            message = <div className="alert alert-danger">
                <p>{this.state.error}</p>
            </div>
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
            </Grid>)
    },

    handleChange: function (field, e) {
        var nextState = {};
        nextState[field] = e.target.checked
        this.setState(nextState)
    },

    handleSubmit: function () {
        var showSuccess = () => {this.transitionTo('/community', {registrationSuccess:true})};
        var showFailure = () => {this.setState({error: "An error occured"})};

        if (this.refs.contactForm.isValid()) {
            var formData = this.refs.contactForm.getFormData();
            this.setState({submitted: formData});
            var url = config.backend_url + '/accounts/update/';

            var fd = new FormData();
            $.each(formData, function (k, v) {
                fd.append(k, v);
            });
            var xhr = new XMLHttpRequest();
            xhr.onload = function () {
                var responseData = JSON.parse(this.response);

                if (this.status == 200 && responseData.success === true) {
                    showSuccess();
                } else {
                    showFailure()
                }
            };
            xhr.open('post', url);
            xhr.setRequestHeader('Authorization', 'JWT ' + auth.token())
            xhr.send(fd);
        } else {
            this.setState({error: "Some information was missing"});
        }
    }
});

var EditProfileForm = React.createClass({
    mixins: [Navigation],
    componentDidMount: function() {
        this.retrieveProfile();
    },
    retrieveProfile: function () {
        var url = config.backend_url + '/accounts/get/';
        var token = auth.token();
        var tokenValue = 'JWT ' + token;
        var saveProfileData = (data) => {
            var imagePreviewUrl = '';
            if (data.user.has_image) {
                imagePreviewUrl = config.backend_url + '/site_media/media/' + data.user['id']
            }
            this.setState({data : data.user, imagePreviewUrl: imagePreviewUrl});
        };
        $.ajax({
            type: 'GET',
            headers: {'Authorization': tokenValue},
            url: url,
            success: function (data) {
                saveProfileData(data);
            },
            error: function (xhr, ajaxOptions, thrownError) {
                this.transitionTo('/signin', {}, {target: '/profile'});
            }.bind(this)
        });
    },
    getInitialState: function () {
        return {errors: {}, data:{}}
    },
    isValid: function () {
        var errors = {};
        if (this.refs.password.getDOMNode().value != this.refs.password_confirm.getDOMNode().value) {
            errors["password_confirm"] = "The passwords don't match"
        }
        this.setState({errors: errors});
        var isValid = true;
        for (var error in errors) {
            isValid = false;
            break;
        }

        return isValid
    },
    getFormData: function () {
        var title = this.refs.titlemd.getDOMNode().checked && this.refs.titlemd.getDOMNode().value ||
            this.refs.titlephd.getDOMNode().checked && this.refs.titlephd.getDOMNode().value ||
            this.refs.titleother.getDOMNode().checked && this.refs.titlecustom.getDOMNode().value;

        var data = {
            image: this.state.file
            , deleteImage : this.state.imageDelete
            , password: this.refs.password.getDOMNode().value
            , password_confirm: this.refs.password_confirm.getDOMNode().value
            , firstName: this.refs.firstName.getDOMNode().value
            , lastName: this.refs.lastName.getDOMNode().value
            , title: title
            , affiliation: this.refs.affiliation.getDOMNode().value
            , institution: this.refs.institution.getDOMNode().value
            , city: this.refs.city.getDOMNode().value
            , state: this.refs.state.getDOMNode().value
            , country: this.refs.country.getDOMNode().value
            , phone_number: this.refs.phone_number.getDOMNode().value
            , comment: this.refs.comment.getDOMNode().value
            , include_me: this.refs.include_me.getDOMNode().checked
            , email_me: this.refs.email_me.getDOMNode().checked
            , hide_number: this.refs.hide_number.getDOMNode().checked
            , hide_email: this.refs.hide_email.getDOMNode().checked
        };
        return data
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
        reader.readAsDataURL(file)
    },
    render: function () {
    return <div className="form-horizontal">
            {this.renderImageUpload('image', 'Profile picture')}
            {this.renderPassword('password', 'Password')}
            {this.renderPassword('password_confirm', 'Confirm Password')}
            {this.renderTextInput('firstName', 'First Name', this.state.data.firstName)}
            {this.renderTextInput('lastName', 'Last Name', this.state.data.lastName)}
            {this.renderRadioInlines('title', '', {
                values: [{name: 'M.D.', ref: 'md'}, {name: 'Ph.D', ref: 'phd'}, {name: 'Other', ref: 'other'}]
                , defaultCheckedValue: this.state.data.title
            })}
            {this.renderSelect('affiliation', 'Affiliation', AFFILIATION, this.state.data.affiliation)}
            {this.renderTextInput('institution', 'Institution, Hospital or Company')}
            {this.renderTextInput('city', 'City', this.state.data.city)}
            {this.renderTextInput('state', 'State or Province', this.state.data.state)}
            {this.renderTextInput('country', 'Country', this.state.data.country)}
            {this.renderTextInput('phone_number', 'Phone number', this.state.data.phone_number)}
            {this.renderTextarea('comment', 'Comment', this.state.data.comment)}
            {this.renderCheckBox('include_me', "Include me in the community page", this.state.data.include_me)}
            {this.renderCheckBox('email_me', "Include me in the mailing list", this.state.data.email_me)}
            {this.renderCheckBox('hide_number', "Don't display my phone number on this website", this.state.data.hide_number)}
            {this.renderCheckBox('hide_email', "Don't display my email on this website",this.state.data.hide_email)}
        </div>
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
                              style={{'maxHeight':'160px', 'maxWidth':'160px'}}/></div>
                    <div ><Button bsStyle="link" onClick={handleImageDelete}>Remove picture</Button></div>
                </div>
            );
        }
        if (imageTooBig) {
            error = <p className="bg-danger">Please choose an image less than 4MB</p>
        }
        return this.renderField(id, label,
            <div>
                <input onChange={this.handleImageChange} type="file" accept="image/*"/>
                {imagePreview}
                {error}
            </div>)
    },
    renderTextInput: function (id, label, defaultValue) {
        var handleChange = () => {var oldData = this.state.data; oldData[id]=this.refs[id].value; this.setState({data: oldData})};
        return this.renderField(id, label,
            <input type="text" className="form-control" id={id} ref={id} value={defaultValue} onChange={handleChange}/>
        )
    },
    renderPassword: function (id, label) {
        return this.renderField(id, label,
            <input type="password" className="form-control" id={id} ref={id}/>
        )
    },
    renderTextarea: function (id, label,defaultValue) {
        var handleChange = () => {var oldData = this.state.data; oldData[id]=this.refs[id].value; this.setState({data: oldData})};
        return this.renderField(id, label,
            <textarea className="form-control" id={id} ref={id} value={defaultValue} onChange={handleChange}/>
        )
    },
    renderSelect: function (id, label, values, defaultValue) {
        var options = values.map(function (value) {
            var selected=defaultValue===value;
            return <option key={id+value} value={value} selected={selected}>{value}</option>
        });
        return this.renderField(id, label,
            <select className="form-control" id={id} ref={id}>
                {options}
            </select>
        )
    },
    renderRadioInlines: function (id, label, kwargs) {
        var handleTextChange = () => {var oldData = this.state.data; oldData[id]=this.refs.titlecustom.getDOMNode().value; this.setState({data: oldData})};
        var otherValue = kwargs.defaultCheckedValue;
        var options = kwargs.values.map(function (value) {
            var handleRadioChange = () => {var oldData = this.state.data; oldData[id] = value.name; this.setState({data: oldData})};
            var defaultChecked = false;
            if (value.name == kwargs.defaultCheckedValue) {
                defaultChecked = true;
                otherValue = '';
            }
            if (value.name == 'Other' && otherValue !=='') {defaultChecked=true}
            return <label className="radio-inline">
                <input type="radio" ref={id+value.ref} name={id} value={value.name} checked={defaultChecked} onChange={handleRadioChange}/>
                {value.name}
            </label>;
        }.bind(this));
        options = <span className="col-sm-9">{options}</span>;
        var other =
            <span className="col-sm-3">
            <input className="form-control" type="text" ref="titlecustom" name="titlecustom" value={otherValue} onChange={handleTextChange}/>
            </span>;
        var optionsWithOther = {options, other};
        return this.renderField(id, label, optionsWithOther)
    },
    renderCheckBox: function (id, label, defaultValue) {
        var handleChange = () => {var oldData = this.state.data; oldData[id]=this.refs[id].checked; this.setState({data: oldData})};
        var checkbox = (<label className="radio-inline">
            <input type='checkbox' ref={id} checked={defaultValue} onChange={handleChange} />
            {label}
        </label>);
        return this.renderField(id, "", checkbox);
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
    Profile: Profile
});
