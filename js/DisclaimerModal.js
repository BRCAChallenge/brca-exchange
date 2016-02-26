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
    render() {
        var modalTrigger = React.cloneElement(this.props.children,{onClick:this.open});
        return (
            <div style={{display: "inline", cursor: "pointer"}}>
                {modalTrigger}
                {this.state.showModal ?
                    <Modal onHide={this.close}>
                        <RawHTML html={content.pages.disclaimer} />
                        <div className = "close-button">
                            <Button onClick={this.close}>close</Button>
                        </div>
                    </Modal> : null }
            </div>
        );
    }
});

module.exports = DisclaimerModal;