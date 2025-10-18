/*global require: false, module: false, __dirname: false */
'use strict';
var HtmlWebpackPlugin = require('html-webpack-plugin');
var webpack = require('webpack');
var path = require('path');
var mdih = require('./loaders/markdown-id-headers');

var databaseKey = require('./databaseKey');
var keyParam = databaseKey.map(function (k) { return "key[]=" + k; }).join("&");
var fs = require('fs');

var port = process.env.BRCAPORT || 8080;

module.exports = {
	historyApiFallback: true,
	entry: "./js/index",
	output: {
		path: "build",
		publicPath: "/",
		filename: "[name].js"
	},
	devServer: {
		port: port
	},
	module: {
		preLoaders: [
		  {
		    test: /\.jsx?$/,
		    loaders: ['eslint'],
		  }
		],
		loaders: [
			{ test: /rx-dom/, loader: "imports?define=>false" },
			{ test: /\.css$/, loader: "style!css" },
			{ test: /\.js$/, exclude: /node_modules/, loader: 'babel-loader?optional=es7.objectRestSpread,optional=runtime,cacheDirectory=true'},
			{
				test: /\.(jpe?g|png|gif|svg|eot|woff2?|ttf)(\?v=[0-9]\.[0-9]\.[0-9])?$/i,
				loaders: ['url?limit=10000'],
				exclude: [path.resolve(__dirname, "js/img/favicon")]
			},
			{ test: /\.md/, type: 'asset/source' },
			// This is a custom loader for the database tsv file that emits a compact
			// json file (no repeated property names), and does a simple sanity check,
			// ensuring that the primary key is, in fact, unique.
			{ test: /enigma-database.tsv$/, loader: "url?limit=10000&name=[hash].json!dsv?" + keyParam}
		]
	},
	plugins: [
		new HtmlWebpackPlugin({
			title: "Template Project",
			filename: "index.html",
			template: "page.template"
		}),
		new webpack.OldWatchingPlugin()
	],
	resolve: {
		alias: {
			rx$: 'rx/dist/rx',
			'rx.binding$': 'rx/dist/rx.binding',
			'rx.async$': 'rx/dist/rx.async',
			'rx.coincidence$': 'rx/dist/rx.coincidence',
			// The npm version of d3-tip is broken. We over-ride the transitive
			// dependency here, with the version we've installed.
			'd3-tip$': require.resolve('d3-tip'),

            // FAISAL: needed for isotope, ref: http://isotope.metafizzy.co/extras.html#webpack
            'masonry': 'masonry-layout',
            'isotope': 'isotope-layout'
		},
		extensions: ['', '.js', '.json', '.coffee'],
		root: __dirname + "/js"
	},
	resolveLoader: {
		alias: {
			'copy': 'file-loader?name=[path][name].[ext]&context=./js'
		}
	},
	'markdown-it': {
		html: true,
		use: [mdih]
	}
};
