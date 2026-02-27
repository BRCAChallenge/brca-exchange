/*global module: false, require: false */
'use strict';

import React from 'react';
import { Modal, Button } from 'react-bootstrap';
import RawHTML from './RawHTML';
import content from './content';

class DisclaimerModal extends React.Component {
    state = { showModal: false };
    
    close = () => {
        this.setState({ showModal: false });
    };

    open = () => {
        this.setState({ showModal: true });
    };

    disableResearchMode = () => {
        if(this.props.onToggleMode) this.props.onToggleMode();
        localStorage.setItem('research-mode', false);
    };

    agree = () => {
        if (this.props.onToggleMode) this.props.onToggleMode();
        localStorage.setItem('research-mode', true);
        this.close();
    };

    buttonModal = () => {
        if(localStorage.getItem('research-mode') === 'true') {
            return (
                <span>
		    <Button className="btn-default" onClick={this.disableResearchMode}>Show Summary View of this Variant</Button>
                </span>
            );
        } else {
            return (
                <span>
                    <Button className="btn-default" onClick={this.open}>{this.props.text}</Button>

                    {
                        this.state.showModal &&
                        <Modal show={true} onHide={this.close}>
			    <RawHTML html={content.pages.researchWarning}/>
			    <div className="d-flex gap-2">
                                <Button onClick={this.agree}>YES</Button>
                                <Button onClick={this.close}>NO</Button>
			    </div>
                        </Modal>
                    }
                </span>
            );
        }
    };

    linkModal = () => {
        return (
            <span>
                <a href="#" onClick={(e) => { e.preventDefault(); this.open(); }}>{this.props.text}</a>
		{this.state.showModal && (
                    <Modal show={true} onHide={this.close}>
                        <RawHTML html={content.pages.disclaimer} />
			<div className="d-flex gap-2">
                            <Button onClick={this.close}>OK</Button>
			</div>
                    </Modal>
		)}
            </span>
        );
    };

    render() {
        if(this.props.buttonModal) {
            return this.buttonModal();
        } else {
            return this.linkModal();
        }
    }
}

export default DisclaimerModal;
