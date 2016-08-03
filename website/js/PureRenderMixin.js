/*global require: false, module: false */
'use strict';

var _ = require('underscore');

// deep equality check for pure rendering
module.exports = {
  shouldComponentUpdate: function(nextProps, nextState) {
    return !_.isEqual(this.props, nextProps) ||
           !_.isEqual(this.state, nextState);
  }
};
