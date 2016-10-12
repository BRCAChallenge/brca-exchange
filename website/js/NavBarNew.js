'use strict';

var React = require('react');
var config = require('./config');
var classNames = require('classnames');


var {Navbar, Nav, DropdownButton} = require('react-bootstrap');
var {Link} = require('react-router');

var NavLink = React.createClass({
    render: function () {
        var {children, ...otherProps} = this.props;
        return (
            <li>
                <Link {...otherProps} role='button'>
                    {children}
                </Link>
            </li>
        );
    }
});

var NavBarNew = React.createClass({
    close: function () {
        this.refs.about.setState({open: false});
    },
    getInitialState: function () {
        return {
            showModal: false,
            isBeta: config.environment === "beta"
        };
    },
    shouldComponentUpdate: function (nextProps, nextState) {
        // Only rerender if path has change or the research mode changes, ignoring query.
        var d3TipDiv = document.getElementsByClassName('d3-tip-selection');
        if (d3TipDiv.length !== 0 && d3TipDiv[0].style.opacity !== '0') {
            d3TipDiv[0].style.opacity = '0';
            d3TipDiv[0].style.pointerEvents = 'none';
        }
        return this.props.mode !== nextProps.mode ||
            this.state.loggedin !== nextState.loggedin ||
            this.props.path.split(/\?/)[0] !== nextProps.path.split(/\?/)[0];
    },
    activePath: function (path, tab) {
        var navPath = (path === "") ? "home" : path.split("/")[0];
        return ((navPath === tab) ? "active" : "");
    },
    render: function () {
        var {path} = this.props;
        var brand = (
            <a className="navbar-brand" href="/">
                <h1>
                    <span className="BRCA">BRCA</span>
                    <span className="exchange"> Exchange</span>
                </h1>
                {this.props.mode === 'research_mode' && <span id="research-label" className="label label-info">All Public Data</span>}
                {this.props.mode === 'default' && <span id="research-label" className="label label-info">Expert Reviewed</span>}
            </a>);
        return (
            <div className={classNames("navbar-container", {"beta": this.state.isBeta})}>
                <Navbar fixedTop brand={brand} toggleNavKey={0}>
                    <Nav eventKey={0} navbar right>
                        <NavLink to='/'>Home</NavLink>
                        <NavLink to='/variants'>Variants</NavLink>
                        <NavLink to='/community'>Community</NavLink>
                        <NavLink to='/help'>Help</NavLink>
                        <DropdownButton className={this.activePath(path, "about")} ref='about' title='More'>
                            <NavLink onClick={this.close} to='/about/variation'>
                                BRCA1, BRCA2, and Cancer
                            </NavLink>
                            <NavLink onClick={this.close} to='/about/history'>
                                History of the BRCA Exchange
                            </NavLink>
                            <NavLink onClick={this.close} to='/about/thisSite'>
                                This Site
                            </NavLink>
                            <NavLink onClick={this.close} to='/releases'>
                                Previous Data Releases
                            </NavLink>
                            <NavLink onClick={this.close} to='/about/api'>
                                Webservices for API Data Access
                            </NavLink>
                        </DropdownButton>
                    </Nav>
                    {this.state.isBeta && <div className='beta-header'>This is a beta version of the BRCA Exchange. Please note that some variant information and website features displayed here are under review - for the most up-to-date finalized information, please refer to <a href="http://brcaexchange.org">www.brcaexchange.org</a>. If you encounter any issues while using the beta website, please report them to <a href="mailto:brcaexchange@gmail.com">brcaexchange@gmail.com</a>.</div>}
                </Navbar>
            </div>
        );
    }
});


module.exports = {
    NavLink,
    NavBarNew
};
