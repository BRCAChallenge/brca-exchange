/*global module: false */
'use strict';

var slugify = require('../js/slugify');

module.exports = function(md) {
    var originalHeadingOpen = md.renderer.rules.heading_open;
    // looks for a pattern like this: ((more than one of any char that's not a paren))
    var explicitIDRegex = /\(\(([^)]+)\)\)/;

    md.renderer.rules.heading_open = function (tokens, idx, options, env, self) {
        tokens[idx].attrs = tokens[idx].attrs || [];

        var title = tokens[idx + 1].children.reduce(function (acc, t) {
            return acc + t.content;
        }, '');

        const foundExplictID = explicitIDRegex.exec(title);
        const chosenID = slugify(foundExplictID ? foundExplictID[1] : title);

        // clean up the explicit ID text, if we used it
        if (foundExplictID) {
            // and remove it from the text itself
            tokens[idx + 1].children.map(function (c) {
                c.content = c.content.replace(explicitIDRegex, "");
            })
        }

        tokens[idx].attrs.push(['id', chosenID]);

        if (originalHeadingOpen) {
            return originalHeadingOpen.apply(this, arguments);
        } else {
            return self.renderToken.apply(self, arguments);
        }
    };
};
