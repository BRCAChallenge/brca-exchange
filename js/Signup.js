'use strict';

var React = require('react');
var backend = require('./backend');
var Rx = require('rx');
require('rx-dom');
var _ = require('underscore');
var $ = require('jquery');

var Cookies = require('js-cookie');


var AFFILIATION = [
    'I lead a testing lab',
    'I am a member of a teting lab',
    'I lead a research lab',
    'I am a member of a research lab',
    'I lead an advocacy group',
    'I work at an advocacy group',
    'I am a genetic councelor'];


var Signup = React.createClass({
    getInitialState: function () {
        return {
            submitted: null,
            success: null
        }
    },
    render: function () {
        var message;
        if (this.state.success === true) {
            message = <div className="alert alert-success">
                <p>Account created successfully</p>
            </div>
        } else if (this.state.success === false) {
            message = <div className="alert alert-danger">
                <p>An error occurred creating this account</p>
            </div>
        }
        return <div>
            <div className="panel panel-default">
                {message}
                <div className="panel-heading clearfix">
                    <h3 className="panel-title pull-left">Signup</h3>

                </div>
                <div className="panel-body">
                    <SignupForm ref="contactForm"
                    />
                </div>
                <div className="panel-footer">
                    <button type="button" className="btn btn-primary btn-block" onClick={this.handleSubmit}>Submit
                    </button>
                </div>
            </div>
        </div>
    },

    handleChange: function (field, e) {
        var nextState = {};
        nextState[field] = e.target.checked
        this.setState(nextState)
    },

    handleSubmit: function () {
        if (this.refs.contactForm.isValid()) {
            var formData = this.refs.contactForm.getFormData();

            this.setState({submitted: formData})

            var url = backend.databaseUrl + '/accounts/register/';

            $.ajax({
                url: url,
                data: formData,
                dataType: 'json',
                crossDomain: true,
                method: 'POST',
                success: function (data) {
                    this.setState({success: data.success})
                }.bind(this),
                error: function (xhr, status, err) {
                    this.setState({success: false})

                }.bind(this)
            });
        }
    }
});

var SignupForm = React.createClass({
    getInitialState: function () {
        return {errors: {}}
    },
    isValid: function () {
        var compulsory_fields = ['email', 'email_confirm', 'password', 'password_confirm', 'firstName',
            'lastName', 'city', 'state', 'country', 'phoneNumber'];
        var errors = {};
        if (this.refs.email.getDOMNode().value != this.refs.email_confirm.getDOMNode().value) {
            errors["email_confirm"] = "The emails don't match"
        }
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
        var title = this.refs.titlemd.getDOMNode().checked && this.refs.titlemd.getDOMNode().value ||
            this.refs.titlephd.getDOMNode().checked && this.refs.titlephd.getDOMNode().value ||
            this.refs.titleother.getDOMNode().checked && this.refs.titlecustom.getDOMNode().value;

        var data = {
            email: this.refs.email.getDOMNode().value
            , email_confirm: this.refs.email_confirm.getDOMNode().value
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
            , phoneNumber: this.refs.phoneNumber.getDOMNode().value
            , hideNumber: this.refs.hideNumber.getDOMNode().checked
            , hideEmail: this.refs.hideEmail.getDOMNode().checked
        };
        return data
    },
    render: function () {
        return <div className="form-horizontal">
            {this.renderTextInput('email', 'Email *')}
            {this.renderTextInput('email_confirm', 'Confirm Email *')}
            {this.renderPassword('password', 'Password *')}
            {this.renderPassword('password_confirm', 'Confirm Password *')}
            {this.renderTextInput('firstName', 'First Name *')}
            {this.renderTextInput('lastName', 'Last Name *')}
            {this.renderRadioInlines('title', '', {
                values: [{name: 'M.D.', ref: 'md'}, {name: 'Ph.D', ref: 'phd'}, {name: 'Other', ref: 'other'}]
                , defaultCheckedValue: 'M.D.'
            })}
            {this.renderSelect('affiliation', 'Affiliation', AFFILIATION)}
            {this.renderTextInput('institution', 'Institution, Hospital or Company')}
            {this.renderTextInput('city', 'City *')}
            {this.renderTextInput('state', 'State or Province *')}
            {this.renderTextInput('country', 'Country *')}
            {this.renderTextInput('phoneNumber', 'Phone number *')}
            {this.renderCheckBox('hideNumber', "Don't display my phone number on this website")}
            {this.renderCheckBox('hideEmail', "Don't display my email on this website")}


        </div>
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
    renderTextarea: function (id, label) {
        return this.renderField(id, label,
            <textarea className="form-control" id={id} ref={id}/>
        )
    },
    renderSelect: function (id, label, values) {
        var options = values.map(function (value) {
            return <option value={value}>{value}</option>
        });
        return this.renderField(id, label,
            <select className="form-control" id={id} ref={id}>
                {options}
            </select>
        )
    },
    renderRadioInlines: function (id, label, kwargs) {
        var options = kwargs.values.map(function (value) {
            var defaultChecked = (value.name == kwargs.defaultCheckedValue)
            return <label className="radio-inline">
                <input type="radio" ref={id+value.ref} name={id} value={value.name} defaultChecked={defaultChecked}/>
                {value.name}
            </label>;
        });
        var other = (<label className="radio-inline">
            <input type="text" ref="titlecustom" name="titlecustom"/>
        </label>);
        var optionsWithOther = {options, other};
        return this.renderField(id, label, optionsWithOther)
    },
    renderCheckBox: function (id, label) {
        var checkbox = (<label className="radio-inline">
            <input type='checkbox' ref={id}/>
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
    Signup: Signup
});