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
    agree() {
        localStorage.setItem('research-mode',true);
        this.close();
    },
    render() {
        var modalTrigger = React.cloneElement(this.props.children,{onClick:this.open});
        return (
            <span>
                {modalTrigger}
                {this.state.showModal ?
                    <Modal onRequestHide={this.close}>
                        <RawHTML html={content.pages.disclaimer} />
                            <Button onClick={this.agree}>OK</Button>
                            <Button onClick={this.close}>Cancel</Button>
                    </Modal> : null }
            </span>
        );
    }
});

module.exports = DisclaimerModal;