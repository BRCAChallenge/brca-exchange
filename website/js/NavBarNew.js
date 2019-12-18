'use strict';

var React = require('react');
var config = require('./config');
var classNames = require('classnames');
var content = require('./content');

var RawHTML = require('./RawHTML');
var {Navbar, Nav, DropdownButton, Modal, Button, OverlayTrigger, Popover} = require('react-bootstrap');
var {Link} = require('react-router');


var brcaHeaderLogo = require('./img/brca-logo-transp.png');

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

const ModeButton = React.createClass({
    render: function() {
        const {mode, toggleMode} = this.props;

        const popper = (mode === 'research_mode')
        ? (
            <Popover title="Change View">
            The BRCA Exchange Detail View shows information drawn from multiple databases and is intended to provide professional users a set of annotations which is as comprehensive as possible. For summary information with expert interpretations, click this button to switch to the Summary View.
            </Popover>
        )
        : (
            <Popover title="Change View">
            The BRCA Exchange Summary View shows the clinical significance as reviewed by the expert ENIGMA consortium. For additional variant information, click this button to switch to the Detail View.
            </Popover>
        );

        return (
            <OverlayTrigger placement='bottom' delayShow={300} overlay={popper}>
                <span id="research-label" className="label label-info" style={{cursor: 'help'}} onClick={toggleMode}>
                    {`${mode === 'research_mode' ? "Detail" : "Summary"} View`}
                </span>
            </OverlayTrigger>
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
    toggleMode: function (e) {
        e.preventDefault();
        if (this.props.mode === "research_mode") {
            this.setState({ showModal: false }, function() {
                this.props.toggleMode();
                this.forceUpdate();
            });
        } else {
            this.setState({showModal: true}, function() {
                this.forceUpdate();
            });
        }
    },
    getModal: function () {
        return (
            <Modal onRequestHide={() => this.setState({ showModal: false }, function() {this.forceUpdate();})}>
                <RawHTML html={content.pages.researchWarning}/>
                <Button onClick={() => {this.setState({ showModal: false }, function() {
                    this.props.toggleMode();
                    this.forceUpdate();
                });}}>Yes</Button>
                <Button onClick={() => this.setState({ showModal: false }, function() {this.forceUpdate();})}>No</Button>
            </Modal>
        );
    },
    render: function () {
        const {path} = this.props;
        const brand = (
            <a href="/" className="branding-clickable">
                <img alt="BRCA Exchange Logo" className="logo-img" src={brcaHeaderLogo} height="40" />

                <div className="brand-collapser">
                    <h1>
                        <span className="BRCA">BRCA</span>
                        <span className="exchange"> Exchange</span>
                    </h1>

                    <ModeButton mode={this.props.mode} toggleMode={this.toggleMode} />
                </div>
            </a>
        );

        return (
            <div className={classNames("navbar-container", {"beta": this.state.isBeta})}>
                <Navbar fixedTop brand={brand} toggleNavKey={0}>
                    <Nav eventKey={0} navbar right>
                        <NavLink to='/'>Home</NavLink>
                        <NavLink to='/variants'>Variants</NavLink>
                        <NavLink to='/community'>Community</NavLink>
                        <NavLink to='/help'>Help</NavLink>
                        <DropdownButton className={this.activePath(path, "about")} ref='about' title='More'>
                            <NavLink onClick={this.close} to='/about/thisSite'>
                                This Site
                            </NavLink>
                            <NavLink onClick={this.close} to='/factsheet'>
                                Facts & Stats
                            </NavLink>
                            <NavLink onClick={this.close} to='/releases'>
                                Previous Data Releases
                            </NavLink>
                            <NavLink onClick={this.close} to='/about/api'>
                                Webservices for API Data Access
                            </NavLink>
                            <NavLink onClick={this.close} to='/about/app'>
                                The Mobile App
                            </NavLink>
                            <NavLink onClick={this.close} to='/about/dataSubmissionPolicy'>
                                Data Submission Policy
                            </NavLink>
                            <NavLink onClick={this.close} to='/whydonate'>
                                Donate
                            </NavLink>
                        </DropdownButton>
                    </Nav>
                    {this.state.isBeta && false && <div className='beta-header'>This is a beta version of the BRCA Exchange. Please note that some variant information and website features displayed here are under review - for the most up-to-date finalized information, and to join our community, please refer to <a href="https://brcaexchange.org">www.brcaexchange.org</a>. If you encounter any issues while using the beta website, please report them to <a href="mailto:brcaexchange@gmail.com">brcaexchange@gmail.com</a>.</div>}
                </Navbar>
                {this.state.showModal && this.getModal()}
            </div>
        );
    }
});


module.exports = {
    NavLink,
    NavBarNew
};
