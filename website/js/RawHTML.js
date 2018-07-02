/*global module: false, require: false */
'use strict';

var React = require('react');
var jQuery = require('jquery');

var RawHTML = React.createClass({
    componentDidMount() {
        const {scrollToFragment} = this.props;

        // make all the anchor tags use scrollToFragment, if it was specified
        if (scrollToFragment) {
            jQuery(".markdown a").click(function() {
                const fragment = this.getAttribute("href").slice(1);
                scrollToFragment(fragment);
                history.replaceState(undefined, undefined, '#' + fragment);
            });
        }
    },
    render: function() {
        var {html, ...otherProps} = this.props;
        return (
            <div className='markdown' {...otherProps} dangerouslySetInnerHTML={{__html: html}} />
        );
    }
});

module.exports = RawHTML;
