'use strict';

var React = require('react');
var content = require('./content');
var RawHTML = require('./RawHTML');
var $ = require('jquery');
var config  = require('./config')
var {Grid, Row, Col, Button} = require('react-bootstrap');
var {Navigation, Link} = require('react-router');

var MailingList = React.createClass({
    mixins: [Navigation],
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
        }
        if (this.state.success) {
            return <Grid id="main-grid">
                <Row>
                    <Col sm={10} smOffset={1} md={8} mdOffset={2} className="text-center">
                        <h2>Thank you for signing up for our mailing list.</h2>
                        <h2>You will receive an email asking you to confirm your email address shortly.</h2>
                    </Col>
                </Row>
                <Row>
                    <Col sm={10} smOffset={1} md={8} mdOffset={2} className="text-center">
                        <h4>Would you also like to be listed as part of the BRCA community? This is a great way to show your support as well as make it easier for collaborators and colleagues to find you.</h4>
                    </Col>
                </Row>
                <Row>
                    <Col sm={10} smOffset={1} md={8} mdOffset={2} className="text-center">
                        <Link to="/signup"><Button>Join our Community Space!</Button></Link>
                    </Col>
                </Row>
            </Grid>
        } else {
            return <Grid id="main-grid">
                <Row id="message">
                    {message}
                </Row>
                <Row>
                    <Col sm={10} smOffset={1} className="text-center">
                        <h4>Join the BRCA Exchange mailing list by submitting your email below.</h4>
                        <h4>See <a href="#">About</a> for more information about this initiative.</h4>
                        <br />
                    </Col>
                </Row>
                <Row id="form">
                    <Col md={8} mdOffset={2}>
                        <MailingListForm ref="contactForm"/>
                    </Col>
                </Row>
                <Row id="submit">
                    <Col md={6} mdOffset={3}>
                        <Button type="button" className="btn btn-primary btn-block" onClick={this.handleSubmit}>
                            Sign up for mailing list
                        </Button>
                    </Col>
                </Row>
            </Grid>
        }
    },

    handleChange: function (field, e) {
        var nextState = {};
        nextState[field] = e.target.checked
        this.setState(nextState)
    },

    handleSubmit: function () {
        var showSuccess = () => this.setState({success: true});
        var showFailure = (msg => {this.setState({error: msg})});

        if (this.refs.contactForm.isValid()) {
            var formData = this.refs.contactForm.getFormData();
            this.setState({submitted: formData});
            var url = config.backend_url + '/accounts/mailinglist/';

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
                    var message = responseData.error;
                    if (message === null) {
                        message = "Could not complete mailing list signup";
                    }
                    showFailure(message)
                }
            };
            xhr.open('post', url);
            xhr.send(fd);
        } else {
            this.setState({error: "Some information was missing"});
        }
    }
});

var MailingListForm = React.createClass({
    getInitialState: function () {
        return {errors: {}, file: '', imagePreviewUrl: null, captcha: ""}
    },
    componentDidMount: function() {
        var me = this;
        onRecaptchaLoad(function() {
            grecaptcha.render(me.refs.mailinglistCAPTCHA.getDOMNode(), {sitekey: config.captcha_key, callback: function(resp) {
                me.setState({captcha: resp});
            }});
        });
    },
    isValid: function () {
        var compulsory_fields = ['email', 'email_confirm'];
        var errors = {};
        if (this.refs.email.getDOMNode().value != this.refs.email_confirm.getDOMNode().value) {
            errors["email_confirm"] = "The emails don't match"
        }
        if (this.state.captcha == "") {
            errors["captcha"] = "No CAPTCHA entered"
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
            email: this.refs.email.getDOMNode().value
            , email_confirm: this.refs.email_confirm.getDOMNode().value
            , captcha: this.state.captcha
        };
        return data
    },
    render: function () {
        return <div className="form-horizontal">
            {this.renderTextInput('email', 'Email *')}
            {this.renderTextInput('email_confirm', 'Confirm Email *')}
            {this.renderCAPTCHA('captcha','CAPTCHA *')}
        </div>
    },
    renderTextInput: function (id, label) {
        return this.renderField(id, label,
            <input type="text" className="form-control" id={id} ref={id}/>
        )
    },
    renderCAPTCHA: function(id, label) {
        return this.renderField(id, label, <div ref="mailinglistCAPTCHA"></div>);
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

var trim = function () {
    var TRIM_RE = /^\s+|\s+$/g
    return function trim(string) {
        return string.replace(TRIM_RE, '')
    }
}();

function $c(staticClassName, conditionalClassNames) {
    var classNames = []
    if (typeof conditionalClassNames == 'undefined') {
        conditionalClassNames = staticClassName
    }
    else {
        classNames.push(staticClassName)
    }
    for (var className in conditionalClassNames) {
        if (!!conditionalClassNames[className]) {
            classNames.push(className)
        }
    }
    return classNames.join(' ')
}

module.exports = ({
    MailingList: MailingList,
    trim:  trim,
    $c : $c
});
