var webpack = require('webpack');
var path = require('path');
var config = require('./webpack.config');

config.output.filename = "testBundle.js";
config.output.path = path.resolve(__dirname, 'build');
config.entry = './test/all.js';
module.exports = config;
