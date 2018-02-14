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

var {Grid, Panel} = require('react-bootstrap');
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
                <Panel
                    header={faq.question}
                    collapsable={true}
                    defaultExpanded={this.state.selectedFAQ === question}
                    onSelect={(event) => this._showFaq(question, event)}>
                    <RawHTML html={faq.content} />
                </Panel>
            );
        }, this);
    },
    _showFaq: function(question, event) {
        // stops page from scrolling to top and clearing the url
        event.preventDefault();

        // sets url for most recently selected question
        this.transitionTo(`/faq#${question}`);

        this.setState({ selectedFAQ: question });
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

