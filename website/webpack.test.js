var webpack = require('webpack');
var config = require('./webpack.config');

config.output.filename = "testBundle.js";
config.output.path = "./build";
config.entry = './test/all.js';
module.exports = config;
