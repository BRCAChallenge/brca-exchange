/*eslint-env browser */
/*global require: false */
'use strict';

var React = require('react');

// shims for older browsers
require('babel/polyfill');
require('es5-shim');
require('es5-shim/es5-sham');

var RawHTML = require('./RawHTML');
require('bootstrap/dist/css/bootstrap.css');

var slugify = require('./slugify');
var content = require('./content');

var {Grid, Col, Row} = require('react-bootstrap');
var {Navigation} = require('react-router');

const navbarHeight = 70; // XXX This value MUST match the setting in custom.css

var Faq = React.createClass({
    mixins: [Navigation],
    getInitialState() {
        return { selectedFAQ: slugify(window.location.hash.slice(1)) };
    },
    _getFaqs: function() {
        return content.faqs.map(function(faq) {
            let question = slugify(faq.question);
            return (
                <Row>
                    <Col smOffset={1} sm={10}>
                        <h3 className="faq-question" onClick={() => this._showFaq(question)}>{faq.question}</h3>
                    </Col>
                    <Col key={question} smOffset={1} sm={10} className={this.state.selectedFAQ === question ? 'show' : 'hide'}>
                        <RawHTML html={faq.content} />
                    </Col>
                </Row>
            );
        }, this);
    },
    _showFaq: function(question) {
        let slugifiedQuestion = slugify(question);
        this.transitionTo(`/faq#${slugifiedQuestion}`);
        this.setState({ selectedFAQ: slugifiedQuestion });
    },
    componentDidMount: function () {
        if (this.state.selectedFAQ !== '') {
            let that = this;
            setTimeout(function () {
                var el = document.getElementById(that.state.selectedFAQ);
                if (el) {
                    window.scrollTo(0, el.getBoundingClientRect().top - navbarHeight);
                }
            }, 0);
        }
    },
    render: function () {
        return (
            <Grid id="main-grid" className="faq">
                {this._getFaqs()}
            </Grid>
        );
    }
});

// {selectedFAQ === '' ? null :
//             <style>{`#${selectedFAQ} { animation-name: emphasis; animation-duration: 10s; } `}</style>}

module.exports = Faq;

