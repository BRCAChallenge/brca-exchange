/*global module: false, require: false */
'use strict';

const React = require('react');
const $ = require('jquery');
import {idHelpClicked} from "./util";

const RawHTML = React.createClass({
    componentDidMount() {
        const {scrollToFragment} = this.props;

        // make all the anchor tags use scrollToFragment, if it was specified
        if (scrollToFragment) {
            $(".markdown a").click(function() {
                const fragment = this.getAttribute("href").slice(1);
                scrollToFragment(fragment);
                history.replaceState(undefined, undefined, '#' + fragment);
            });
        }

        // add tooltips to anything that has an id
        const $idElems = $("*[id]", this.refs.me.getDOMNode());
        $idElems.each((idx, elem) => {
            $(elem).addClass("id-associated");
            $("<a>")
                .addClass("id-helper-tip")
                .html("&para;")
                .click((e) => { idHelpClicked(e); })
                .attr("href", `#${$(elem).attr("id")}`)
                .appendTo(elem);
        });
    },
    render: function() {
        var {html, ...otherProps} = this.props;
        return (
            <div ref="me" className='markdown' {...otherProps} dangerouslySetInnerHTML={{__html: html}} />
        );
    }
});

module.exports = RawHTML;
