/*global module: false, require: false */
'use strict';

import {HARDLINK_GLYPH, idHelpClicked} from "./components/help/HardlinkHelper";

const React = require('react');
const ReactDOM = require('react-dom');
const $ = require('jquery');

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

        if (this.props.hardlinks) {
            // add tooltips to anything that has an id
            const $idElems = $("*[id]", ReactDOM.findDOMNode(this.refs.me));
            $idElems.each((idx, elem) => {
                // TODO: ideally we'd replace this ad-hoc html element with the actual HardlinkHelper, although that
                //  would be tricky/hacky since the contents of RawHTML aren't part of the react component tree...
                //  (for now, it resembles the actual version enough that it shouldn't be a problem.)
                $(elem).addClass("id-associated");
                $("<a>")
                    .addClass("id-helper-tip")
                    .html(HARDLINK_GLYPH)
                    .click((e) => { idHelpClicked(e); })
                    .attr("href", `#${$(elem).attr("id")}`)
                    .appendTo(elem);
            });
        }
    },
    render: function() {
        var {html, ...otherProps} = this.props;
        return (
            <div ref="me" className='markdown' {...otherProps} dangerouslySetInnerHTML={{__html: html}} />
        );
    }
});
RawHTML.propTypes = {
    html: React.PropTypes.string.isRequired,
    hardlinks: React.PropTypes.bool
};

module.exports = RawHTML;
