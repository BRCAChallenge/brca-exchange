/*global require: false, module: false */
'use strict';

var webpack = require('webpack');
var config = require('./webpack.config');

config.output.filename = "[name].[contenthash].js";
config.output.chunkFilename = "[contenthash].bundle.js";

config.optimization = {
    ...config.optimization,
    minimize: true,
    moduleIds: 'deterministic',
};

config.plugins = config.plugins.concat([
    new webpack.DefinePlugin({
    	'process.env.NODE_ENV': JSON.stringify('production')
    })
]);

module.exports = config;
