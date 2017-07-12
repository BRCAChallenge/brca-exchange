/*global __dirname: false */
const webpack = require('webpack');
const path = require('path');
const config = require('./webpack.config');

config.output.filename = "testBundle.js";
config.output.path = path.resolve(__dirname, "build");
config.entry = './test/all.js';
module.exports = config;
