/*eslint-env browser */
/*global require: false */
'use strict';

// shims for older browsers
require('babel/polyfill');
require('es5-shim');
require('es5-shim/es5-sham');

// library CSS includes
require('bootstrap/dist/css/bootstrap.css');
require('font-awesome-webpack');
// local CSS/image includes
require('css/bootstrap-xlgrid.css'); // adds xl, xxl, xxxl grid sizes to bootstrap 3
require('css/custom.css');
require('./data/favicons');

// library includes
var React = require('react');
var {State, Route, RouteHandler, HistoryLocation, run, DefaultRoute} = require('react-router');

if (typeof console === "undefined") {
    window.console = {
        log: function () {}
    };
}

// partials
var {NavBarNew} = require('./partials/NavBarNew');
var Footer = require('./partials/Footer');

// routes
var Home = require('./routes/Home');
var About = require('./routes/About');
var FactSheet = require('./routes/FactSheet');
var Help = require('./routes/Help');
var Community = require('./routes/Community');
var {Signup} = require('./routes/Signup');
var {Signin, ResetPassword} = require('./routes/Signin');
var {MailingList} = require('./routes/MailingList');
var {Profile} = require('./routes/Profile');
var {ConfirmEmail} = require('./routes/ConfirmEmail');
var {ChangePassword} = require('./routes/ChangePassword');
var Database = require('./routes/Database'); // pseudo-route
var VariantDetail = require('./routes/VariantDetail');
var {Releases, Release} = require('./routes/Releases.js');


var Application = React.createClass({
    mixins: [State],
    onChildToggleMode: function() {
        this.toggleMode();
    },
    getInitialState: function () {
        return {
            mode: (localStorage.getItem("research-mode") === 'true') ? 'research_mode' : 'default',
        };
    },
    toggleMode: function () {
        if (this.state.mode === 'research_mode') {
            localStorage.setItem('research-mode', false);
            this.setState({mode: 'default'});
        } else {
            localStorage.setItem('research-mode', true);
            this.setState({mode: 'research_mode'});
        }
    },
    render: function () {
        var path = this.getPath().slice(1);
        return (
            <div>
                <NavBarNew path={path} mode={this.state.mode} />
                <RouteHandler toggleMode={this.onChildToggleMode} mode={this.state.mode} />
                <Database
                    mode={this.state.mode}
                    toggleMode={this.onChildToggleMode}
                    show={path.indexOf('variants') === 0} />
                <Footer />
            </div>
        );
    }
});

var routes = (
    <Route handler={Application}>
        <DefaultRoute handler={Home}/>
        <Route path='about/:page' handler={About}/>
        <Route path='factsheet' handler={FactSheet}/>
        <Route path='help' handler={Help}/>
        <Route path='community' handler={Community}/>
        <Route path='signup' handler={Signup}/>
        <Route path='signin' handler={Signin}/>
        <Route path='mailinglist' handler={MailingList}/>
        <Route path='reset_password' handler={ResetPassword}/>
        <Route path='profile' handler={Profile}/>
        <Route path='confirm/:activationCode' handler={ConfirmEmail}/>
        <Route path='reset/:resetToken' handler={ChangePassword}/>
        <Route path='variants' />
        <Route path='variant/:id' handler={VariantDetail}/>
        <Route path='releases' handler={Releases}/>
        <Route path='release/:id' handler={Release}/>
    </Route>
);

var main = document.getElementById('main');

run(routes, HistoryLocation, (Root) => {
  React.render(<Root/>, main);
});
