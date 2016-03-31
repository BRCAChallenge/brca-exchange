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
        this.props.onToggleMode.forceUpdate();
        localStorage.setItem('research-mode',false);
    },
    agree() {
        this.props.onToggleMode.forceUpdate();
        localStorage.setItem('research-mode',true);
        this.close();
    },
    buttonModal() {
        if(localStorage.getItem('research-mode') == 'true') {
            return (
                <div className="form-group">
                <Button className="btn-sm" onClick={this.disable_research_mode}>Return to the default view</Button>
                </div>
            );
        } else {
            return (
                <div className="form-group">
                    <Button className="btn-sm" onClick={this.open}>{this.props.text}</Button>
                    {this.state.showModal ?
                        <Modal onRequestHide={this.close}>
                            <RawHTML html={content.pages.researchWarning} />
                            <Button onClick={this.agree}>YES</Button>
                            <Button onClick={this.close}>NO</Button>
                        </Modal> : null }
                </div>
            );
        }
    },
    linkModal() {
        return (
            <span>
                <a style={{cursor:"pointer"}} onClick={this.open}>{this.props.text}</a>
                {this.state.showModal ?
                    <Modal onRequestHide={this.close}>
                        <RawHTML html={content.pages.disclaimer} />
                        <Button onClick={this.close}>OK</Button>
                    </Modal> : null }
            </span>
        );
    },
    render() {
        if(this.props.buttonModal){
            return this.buttonModal()
        } else {
            return this.linkModal()
        }
    }
});

module.exports = DisclaimerModal;