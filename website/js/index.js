/*eslint-env browser */
/*global require: false */

// shims for older browsers
// require('babel/polyfill');
// require('es5-shim');
// require('es5-shim/es5-sham');

// library CSS includes
import 'bootstrap/dist/css/bootstrap.css';
import 'font-awesome-webpack';

// local CSS/image includes
import './css/bootstrap-xlgrid.css'; // adds xl, xxl, xxxl grid sizes to bootstrap 3
import './css/custom.css';
import './data/favicons';

// library includes
import React from 'react';

import { State, Route, RouteHandler, HistoryLocation, run, DefaultRoute } from 'react-router';

// partials
import { NavBarNew } from './partials/NavBarNew';
import Footer from './partials/Footer';

// routes
import Home from './routes/Home';
import About from './routes/About';
import FactSheet from './routes/FactSheet';
import Help from './routes/Help';
import Community from './routes/Community';
import { Signup } from './routes/Signup';
import { Signin, ResetPassword } from './routes/Signin';
import { MailingList } from './routes/MailingList';
import Profile from './routes/Profile';
import ConfirmEmail from './routes/ConfirmEmail';
import ChangePassword from './routes/ChangePassword';
import Database from './routes/Database'; // pseudo-route
import VariantDetail from './routes/VariantDetail';
import { Releases, Release } from './routes/Releases.js';


// no-op log function when there's no console
if (typeof console === "undefined") {
    window.console = {
        log: function () {}
    };
}


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
