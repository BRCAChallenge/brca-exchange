/*global require: false, module: false, __dirname: false */
'use strict';
var HtmlWebpackPlugin = require('html-webpack-plugin');
var webpack = require('webpack');
var path = require('path');

module.exports = {
	historyApiFallback: true,
	entry: "./js/bogorouter",
	output: {
		path: "build",
		publicPath: "/",
		filename: "[name].js"
	},
	module: {
		loaders: [
			{ test: /rx-dom/, loader: "imports?define=>false" },
			{ test: /\.css$/, loader: "style!css" },
			{ test: /\.js$/, exclude: /node_modules/, loader: 'babel-loader?optional=es7.objectRestSpread,optional=runtime,cacheDirectory=true'},
			{
				test: /\.(jpe?g|png|gif|svg|eot|woff2?|ttf)$/i,
				loaders: ['url?limit=10000'],
				exclude: [path.resolve(__dirname, "js/img/favicon")]
			}
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
			'rx.coincidence$': 'rx/dist/rx.coincidence'
		},
		extensions: ['', '.js', '.json', '.coffee'],
		root: __dirname + "/js"
	},
	resolveLoader: {
		alias: {
			'copy': 'file-loader?name=[path][name].[ext]&context=./js',
		}
	}
};
