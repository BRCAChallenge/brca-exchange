/*global module: false */
'use strict';

var slugify = require('../js/slugify');

module.exports = function(md) {
	var originalHeadingOpen = md.renderer.rules.heading_open;
	md.renderer.rules.heading_open = function (tokens, idx, something, somethingelse, self) {
		tokens[idx].attrs = tokens[idx].attrs || [];

		var title = tokens[idx + 1].children.reduce(function (acc, t) {
			return acc + t.content;
		}, '');

		var slug = slugify(title);
		tokens[idx].attrs.push(['id', slug]);

		if (originalHeadingOpen) {
			return originalHeadingOpen.apply(this, arguments);
		} else {
			return self.renderToken.apply(self, arguments);
		}
	};
};
