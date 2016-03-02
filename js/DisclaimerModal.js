/*global module: false, require: false */
'use strict';

var React = require('react');
var {Modal, Button} = require('react-bootstrap');
var RawHTML = require('./RawHTML');
var content = require('./content');

var DisclaimerModal = React.createClass({
    getInitialState() {
        return { showModal: false };
    },
    close() {
        this.setState({ showModal: false });
    },
    open() {
        this.setState({ showModal: true });
    },
    disable_research_mode() {
        localStorage.setItem('research-mode',false);
    },
    agree() {
        localStorage.setItem('research-mode',true);
        this.close();
    },
    research_mode() {
        if(localStorage.getItem('research-mode') == 'true') {
            return <Button onClick={this.disable_research_mode}>Return to the default view</Button>
        } else {
            return (
                <span>
                    <Button onClick={this.open}>Research information on this variant</Button>
                    {this.state.showModal ?
                        <Modal onRequestHide={this.close}>
                            <RawHTML html={content.pages.disclaimer} />
                            <Button onClick={this.agree}>OK</Button>
                            <Button onClick={this.close}>Cancel</Button>
                        </Modal> : null }
                </span>
            );
        }
    },
    general_mode() {
        return (
            <span>
                <a style={{cursor:"pointer"}} onClick={this.open}>disclaimer</a>
                {this.state.showModal ?
                    <Modal onRequestHide={this.close}>
                        <RawHTML html={content.pages.disclaimer} />
                        <Button onClick={this.agree}>OK</Button>
                    </Modal> : null }
            </span>
        );
    },
    render() {
        if(this.props.research_mode){
            return this.research_mode()
        } else {
            return this.general_mode()
        }
    }
});

module.exports = DisclaimerModal;