'use strict';

var React = require('react');
var {State} = require('react-router');
var {Grid, Col, Row} = require('react-bootstrap');

var content = require('../data/content');
var {navbarHeight} = require('../data/constants');
var slugify = require('../helpers/slugify');
var RawHTML = require("../helpers/RawHTML");

var Help = React.createClass({
    mixins: [State],
    componentDidMount: function () {
        var fragment = slugify(window.location.hash.slice(1));
        if (fragment !== '') {
            setTimeout(function () {
                var el = document.getElementById(fragment);
                if (el) {
                    window.scrollTo(0, el.getBoundingClientRect().top - navbarHeight);
                }
            }, 0);
        }
    },
    render: function () {
        var fragment = slugify(window.location.hash.slice(1));
        var helpContent;
        if (localStorage.getItem("research-mode") === 'true') {
            helpContent = content.pages.helpResearch;
        } else {
            helpContent = content.pages.help;
        }
        return (
            <Grid id="main-grid" className="help">
                {fragment === '' ? null :
                    <style>{`#${fragment} { animation-name: emphasis; animation-duration: 10s; } `}</style>}
                <Row>
                    <Col smOffset={1} sm={10}>
                        <RawHTML ref='content' html={helpContent}/>
                    </Col>
                </Row>
            </Grid>
        );
    }
});

module.exports = Help;
