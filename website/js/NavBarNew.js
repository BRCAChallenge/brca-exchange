'use strict';

import React from 'react';
import config from './config';
import classNames from 'classnames';
import content from './content';

import RawHTML from './RawHTML';
import { Navbar, NavDropdown, Nav, Modal, Button, OverlayTrigger, Popover } from 'react-bootstrap';
import { Link } from 'react-router-dom';


var brcaHeaderLogo = require('./img/brca-logo-transp.png');

class NavLink extends React.PureComponent {
    render() {
        const {children, ...otherProps} = this.props;
	return (
	    <Nav.Link as={Link} {...otherProps} role="button">
		{children}
	    </Nav.Link>
        );
    }
}

class ModeButton extends React.PureComponent {
    render() {
        const {mode, toggleMode} = this.props;

        const popper = (mode === 'research_mode')
        ? (
            <Popover id="change-view-popover">
		<Popover.Header as="h3">Change View</Popover.Header>
		<Popover.Body>
            	    The BRCA Exchange Detail View shows information drawn from multiple databases and is intended to provide professional users a set of annotations which is as comprehensive as possible. For summary information with expert interpretations, click this button to switch to the Summary View.
                </Popover.Body>
	    </Popover>
        )
        : (
            <Popover id="change-view-popover">
		<Popover.Header as="h3">Change View</Popover.Header>
		<Popover.Body>
            The BRCA Exchange Summary View shows the clinical significance as reviewed by the expert ENIGMA consortium. For additional variant information, click this button to switch to the Detail View.
            	</Popover.Body>
	    </Popover>
        );

        return (
            <OverlayTrigger placement='bottom' delay={{ show: 300}} overlay={popper}>
                <span id="research-label" className="label label-info" style={{cursor: 'help'}} onClick={toggleMode}>
                    {`${mode === 'research_mode' ? "Detail" : "Summary"} View`}
                </span>
            </OverlayTrigger>
        );
    }
}

class NavBarNew extends React.Component {
    // no-op: RB v2 dropdowns auto-close on item click
    close = () => {};

    state = {
        showModal: false,
        isBeta: config.environment === "beta"
    };

    shouldComponentUpdate(nextProps, nextState) {
        // Only rerender if path has change or the research mode changes, ignoring query.
        var d3TipDiv = document.getElementsByClassName('d3-tip-selection');
        if (d3TipDiv.length !== 0 && d3TipDiv[0].style.opacity !== '0') {
            d3TipDiv[0].style.opacity = '0';
            d3TipDiv[0].style.pointerEvents = 'none';
        }
        return this.props.mode !== nextProps.mode ||
            this.state.loggedin !== nextState.loggedin ||
            this.state.showModal !== nextState.showModal ||
            this.props.path.split(/\?/)[0] !== nextProps.path.split(/\?/)[0];
    }
    
    activePath(path, tab) {
        var navPath = (path === "") ? "home" : path.split("/")[0];
        return ((navPath === tab) ? "active" : "");
    }

    toggleMode = (e) => {
        e.preventDefault();
        if (this.props.mode === "research_mode") {
            this.setState({ showModal: false }, () => {
                this.props.toggleMode();
                this.forceUpdate();
            });
        } else {
            this.setState({showModal: true}, () => {
                this.forceUpdate();
            });
        }
    }

    render() {
        const {path} = this.props;
        const brand = (
	    <span className="branding-clickable">
                <img alt="BRCA Exchange Logo" className="logo-img" src={brcaHeaderLogo} height="40" />

                <div className="brand-collapser">
                    <h1>
                        <span className="BRCA">BRCA</span>
                        <span className="exchange"> Exchange</span>
                    </h1>

                    <ModeButton mode={this.props.mode} toggleMode={this.toggleMode} />
                </div>
	     </span>
        );

        return (
            <div className={classNames("navbar-container", {"beta": this.state.isBeta})}>
            	<Navbar bg="light" expand="lg" fixed="top">
          		<Navbar.Brand as={Link} to="/">{brand}</Navbar.Brand>
          		<Navbar.Toggle aria-controls="main-navbar" />
          		<Navbar.Collapse id="main-navbar">
            			<Nav className="ms-auto">
              				<Nav.Link as={Link} to="/">Home</Nav.Link>
              				<Nav.Link as={Link} to="/variants">Variants</Nav.Link>
              				<Nav.Link as={Link} to="/community">Community</Nav.Link>
              				<Nav.Link as={Link} to="/help">Help</Nav.Link>
              				<NavDropdown id="about-dropdown" title="More" className={this.activePath(path || '', 'about')}>
                			<NavDropdown.Item as={Link} to="/about/thisSite" onClick={this.close}>This Site</NavDropdown.Item>
                			<NavDropdown.Item as={Link} to="/factsheet" onClick={this.close}>Facts &amp; Stats</NavDropdown.Item>
                			<NavDropdown.Item as={Link} to="/releases" onClick={this.close}>Data Releases</NavDropdown.Item>
                			<NavDropdown.Item as={Link} to="/about/api" onClick={this.close}>Webservices for API Data Access</NavDropdown.Item>
                			<NavDropdown.Item as={Link} to="/about/dataSubmissionPolicy" onClick={this.close}>Data Submission Policy</NavDropdown.Item>
                			<NavDropdown.Item as={Link} to="/whydonate" onClick={this.close}>Donate</NavDropdown.Item>
                			<NavDropdown.Divider />
                			<NavDropdown.Item href="https://brcaexchange.org/blog">Blog</NavDropdown.Item>
              			</NavDropdown>
            			</Nav>
          		</Navbar.Collapse>    
		{this.state.isBeta && false && <div className='beta-header'>This is a beta version of the BRCA Exchange. Please note that some variant information and website features displayed here are under review - for the most up-to-date finalized information, and to join our community, please refer to <a href="https://brcaexchange.org">www.brcaexchange.org</a>. If you encounter any issues while using the beta website, please report them to <a href="mailto:brcaexchange@gmail.com">brcaexchange@gmail.com</a>.</div>}
                </Navbar>

                {
                    this.state.showModal &&
                    <Modal show={true} onHide={() => this.setState({showModal: false})}>
                        <RawHTML html={content.pages.researchWarning}/>
                        <Button onClick={() => {
                            this.setState({showModal: false}, function () {
                                this.props.toggleMode();
                                this.forceUpdate();
                            });
                        }}>Yes</Button>
                        <Button onClick={() => this.setState({showModal: false}, function () {
                            console.log("Modal visible?: " + this.state.showModal);
                        })}>No</Button>
                    </Modal>
                }
            </div>
        );
    }
}

export { NavLink };
export default NavBarNew;
