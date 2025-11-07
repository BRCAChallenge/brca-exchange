/*global require: false, module: false, __dirname: false */
'use strict';

const HtmlWebpackPlugin = require('html-webpack-plugin');
const path = require('path');
// const ESLintPlugin = require('eslint-webpack-plugin'); // optional

const databaseKey = require('./databaseKey');
const keyParam = databaseKey.map(k => `key[]=${k}`).join('&');

const port = process.env.BRCAPORT || 8088;

module.exports = {
  entry: './js/index',
  output: {
    path: path.resolve(__dirname, 'build'),
    publicPath: '/',
    filename: '[name].js',
    clean: false
  },

  devServer: {
    port,
    historyApiFallback: true
  },

  module: {
    rules: [
      // replaces "imports?define=>false"
      {
        test: /rx-dom/,
        use: [{
          loader: 'imports-loader',
          options: { additionalCode: 'var define = false;' }
        }]
      },

      // CSS
      {
        test: /\.css$/,
        use: ['style-loader', 'css-loader']
      },

      // JS/JSX via Babel (now with presets)
      {
        test: /\.(js|jsx)$/,
        exclude: /node_modules/,
        use: {
          loader: 'babel-loader',
          options: {
            cacheDirectory: true,
            presets: [
              ['@babel/preset-env', { modules: false }],
              ['@babel/preset-react', { runtime: 'automatic', development: true }]
            ]
          }
        }
      },

      // Images & fonts
      {
        test: /\.(jpe?g|png|gif|svg|eot|woff2?|ttf)(\?v=\d+\.\d+\.\d+)?$/i,
        type: 'asset',
        parser: { dataUrlCondition: { maxSize: 10 * 1024 } },
        generator: { filename: 'assets/[name][contenthash][ext]' },
        exclude: [path.resolve(__dirname, 'js/img/favicon')]
      },

      // Raw markdown
      { test: /\.md$/, type: 'asset/source' },

      // TSV → emitted compact JSON (keeps your custom pipeline intention)
      {
        test: /enigma-database\.tsv$/,
        type: 'asset/resource',
        generator: { filename: '[contenthash].json' },
        use: [{
          loader: path.resolve(__dirname, 'loaders/tsv-to-json-inline.js')
        }]
      }
    ]
  },

  plugins: [
    new HtmlWebpackPlugin({
      title: 'Template Project',
      filename: 'index.html',
      template: 'page.template'
    })
    // new ESLintPlugin({ extensions: ['js','jsx'], context: path.resolve(__dirname, 'js') })
  ],

  resolve: {
    alias: {
      rx$: 'rx/dist/rx',
      'rx.binding$': 'rx/dist/rx.binding',
      'rx.async$': 'rx/dist/rx.async',
      'rx.coincidence$': 'rx/dist/rx.coincidence',
      'd3-tip$': require.resolve('d3-tip'),
      masonry: 'masonry-layout',
      isotope: 'isotope-layout'
    },
    extensions: ['.js', '.jsx', '.json', '.coffee'],
    modules: [path.resolve(__dirname, 'js'), 'node_modules']
  },

  resolveLoader: {
    alias: {
      copy: 'file-loader?name=[path][name].[ext]&context=./js'
    }
  },

  stats: 'minimal'
};

