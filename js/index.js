/*eslint-env browser */
/*global require: false */
'use strict';

// shims for older browsers
require('babel/polyfill');
require('es5-shim');
require('es5-shim/es5-sham');

require('./favicons');
var React = require('react');
var PureRenderMixin = require('./PureRenderMixin'); // deep-equals version of PRM
require('bootstrap/dist/css/bootstrap.css');
require('font-awesome-webpack');
var Rx = require('rx');
var vcf = require('vcf.js');
require('rx-dom');
require('css/custom.css');
var _ = require('underscore');
//var brcaLogo = require('./img/brca_logo.png');
var logos = require('./logos');
var slugify = require('./slugify');

var content = require('./content');

var databaseUrl = require('../../enigma-database.tsv');
var databaseKey = require('../databaseKey');

var {Grid, Col, Row, Input, Navbar, Nav, Table,
	DropdownButton} = require('react-bootstrap');


var VariantTable = require('./VariantTable');
var VariantSearch = require('./VariantSearch');
var {Navigation, State, Link, Route, RouteHandler,
	HistoryLocation, run, DefaultRoute} = require('react-router');

var navbarHeight = 70; // XXX This value MUST match the setting in custom.css

var variantPathJoin = row => _.map(databaseKey, k => encodeURIComponent(row[k])).join('@@');
var variantPathSplit = id => _.object(databaseKey, _.map(id.split(/@@/), decodeURIComponent));

function readTsv(response) {
	var {header, rows} = JSON.parse(response);
	return {
		records: _.map(rows, row => _.object(header, row))
	};
}

var RawHTML = React.createClass({
	render: function() {
		var {html, ...otherProps} = this.props;
		return (
			<div className='markdown' {...otherProps} dangerouslySetInnerHTML={{__html: html}} />
		);
	}
});

var NavLink = React.createClass({
	render: function () {
		var {children, ...otherProps} = this.props;
		return (
			<li>
				<Link {...otherProps} role='button'>
					{children}
				</Link>
			</li>
		);
	}
});

var NavBarNew = React.createClass({
	close: function () {
		this.refs.about.setState({open: false});
	},
    activePath: function(path, tab) {
        var navPath = (path === "") ? "home" : path.split("/")[0];
        return ((navPath === tab) ? "active" : "");
    },
	render: function () {
        var {path} = this.props;
		var brand = (
			<a className="navbar-brand" href="http://brcaexchange.org">
				<span>
					<b className="BRCA">BRCA</b>
					<span className="exchange"> Exchange</span>
				</span>
			</a>);
		return (
			<div className="navbar-container">
            <Navbar fixedTop brand={brand} toggleNavKey={0}>
				<Nav eventKey={0} navbar right>
					<NavLink to='/'>Home</NavLink>
					<DropdownButton className={this.activePath(path, "about")} ref='about' title='About'>
						<NavLink onClick={this.close} to='/about/history'>
							History of the BRCA Exchange
						</NavLink>
						<NavLink onClick={this.close} to='/about/brca1_2'>
							What are BRCA1 and BRCA2?
						</NavLink>
						<NavLink onClick={this.close} to='/about/variation'>
							BRCA Variation and Cancer
						</NavLink>
					</DropdownButton>
					<NavLink to='/variants'>Variants</NavLink>
					<NavLink to='/help'>Help</NavLink>
				</Nav>
			</Navbar>
            </div>
		);
	}
});

var Footer = React.createClass({
    render: function() {
        return (
            <Grid className="footer">
                <Row><a href="#">disclaimer</a></Row>
                <Row><a href="#">contact us</a></Row>
            </Grid>
        );
    }
});

var Home = React.createClass({
	mixins: [Navigation],
	getInitialState() {
		return {
			index: 0,
			direction: null
		};
	},

	handleSelect(selectedIndex, selectedDirection) {
		this.setState({
			index: selectedIndex,
			direction: selectedDirection
		});
	},
	onSearch(value) {
		this.transitionTo('/variants', null, {search: value});
	},
	render: function() {
		var {suggestions} = this.props;
		var logoItems = _.map(logos, ({id, logo, url}) => (
			<li key={id}><a href={url}>
				<img id={id} src={logo} alt={id + ' logo'} />
			</a></li>
		));
		return (
			<Grid className='home'>
				<Row>
				   	<div className='text-center'>
						<VariantSearch
							id='home-search'
							suggestions={suggestions}
							onSearch={this.onSearch}/>
					</div>
				</Row>
				<Row>
                    <div className="jumbotron">
					    <RawHTML html={content.pages.home} />
				    </div>
                </Row>
				<Row>
			        <div>
                        <ul>
                            <div className="container logo-block">
						        {logoItems}
					        </div>
                        </ul>
                    </div>
				</Row>
			</Grid>
		);
	}
});

var About = React.createClass({
	render: function() {
		var {page} = this.props.params;

		return (
			<Grid>
				<Row>
					<Col md={8} mdOffset={2}>
						<RawHTML html={content.pages[page]} />
					</Col>
				</Row>
			</Grid>
		);
	}
});

var Help = React.createClass({
	mixins: [State],
	componentDidMount: function () {
		var fragment = slugify(window.location.hash.slice(1));
		if (fragment !== '') {
			setTimeout(function () {
				var el = document.getElementById(fragment);
				if (el) {
					window.scrollTo(0, el.getBoundingClientRect().top - navbarHeight);
				}
			}, 0);
		}
	},
	render: function() {
		var fragment = slugify(window.location.hash.slice(1));
		return (
			<Grid className="help">
				{fragment === '' ? null :
					<style>{`#${fragment} { animation-name: emphasis; animation-duration: 10s; } `}</style>}
				<Row>
					<Col md={8} mdOffset={2}>
						<RawHTML ref='content' html={content.pages.help} />
					</Col>
				</Row>
			</Grid>
		);
	}
});

var Database = React.createClass({
	mixins: [Navigation, State, PureRenderMixin],
	showVariant: function (row) {
		this.transitionTo(`/variant/${variantPathJoin(row)}`);
	},
	showHelp: function (title) {
		this.transitionTo(`/help#${slugify(title)}`);
	},
	render: function () {
		var {show, data, suggestions} = this.props,
			{search} = this.getQuery();
		return (
			<Grid style={{display: show ? 'block' : 'none'}}>
				{data ?
					<VariantTable
						ref='table'
						filterValues={{visibleSearch: search || ''}}
						data={data.records}
						suggestions={suggestions}
						keys={databaseKey}
						onHeaderClick={this.showHelp}
						onRowClick={this.showVariant}/>
					: ''}
			</Grid>
		);
	}
});

var MyVariant = React.createClass({ //eslint-disable-line no-unused-vars
	getInitialState: function () {
		return {
			data: null
		};
	},

	render: function() {
		var {data} = this.state;
		var {show} = this.props;
		return (
			<div style={{display: show ? 'block' : 'none'}}>
				<div className="text-center">
					<Input ref='file' type='file' onChange={this.fileChange}/>
				</div>
				<div style={{position: "relative", height: "100px"}}>
					{data ?
						<Row>
							<Col md={10} mdOffset={1}>
								<VariantTable data={data}/>
							</Col>
						</Row>
						: ''}
				</div>
			</div>
		);
	},

	dataReady: function(ev) {
		this.setState({data: vcf.parser()(ev.currentTarget.result)});
	},


	fileChange: function () {
		var file = this.refs.file.getInputDOMNode().files[0],
		reader = new FileReader();
		// XXX This timeout allows the UI to update (close dialog) before loading
		// a potentially large file, which will block the UI.
		// This might also be solved by elminating the animation on Modal close,
		// which is probably the source of the problem.
		window.setTimeout(() => {
			reader.onload = this.dataReady;
			reader.readAsText(file);
		}, 100);
	}
});

var VariantDetail = React.createClass({
	mixins: [Navigation],
	showHelp: function (title) {
		this.transitionTo(`/help#${slugify(title)}`);
	},
	render: function() {
		var {data, params: {id}} = this.props,
			variant = (data && _.findWhere(data.records, variantPathSplit(id))) || {};

		variant = _.omit(variant, ['__HEADER__']);
		var rows = _.map(variant, (v, k) =>
			 <tr key={k}>
				 <td className='help-target'>
					{k}
					<span onClick={() => this.showHelp(k)}
						className='help glyphicon glyphicon-question-sign superscript'/>
				 </td>
				 <td>{v}</td>
			 </tr>);


		return (
			<Grid>
				<Row>
					<div className='text-center Variant-detail-title'>
						<h3>Variant Detail</h3>
					</div>
				</Row>
				<Row>
					<Col md={8} mdOffset={2}>
						<Table striped bordered>
							<tbody>
								{rows}
							</tbody>
						</Table>
					</Col>
				</Row>
			</Grid>
		);
	}
});

var dontSuggest = [
	'Assertion method citation',
	'Citations or URLs for  clinical significance without database identifiers'
];

var flatmap = (coll, fn) => _.flatten(_.map(coll, fn), true);
var minSuggestion = 3;
var rowWords = row => flatmap(_.values(_.omit(row, dontSuggest)), v => v.split(/\s+/));

// Pull out interesting strings from the data, for use in
// auto-completion.
function getSuggestions(data) {
	return _.uniq(flatmap(data, row =>
				_.filter(rowWords(row), w => w.length >= minSuggestion)).sort(),
			true);
}

var Application = React.createClass({
	mixins: [State],
	getInitialState: function () {
		return {data: null};
	},
	componentWillMount: function (){
		Rx.DOM.get(databaseUrl).subscribe(xhr => {
			var data = readTsv(xhr.responseText);
			this.setState({data: data, suggestions: getSuggestions(data.records)});
		});
	},
	render: function () {
		var {data, suggestions} = this.state;
		var path = this.getPath().slice(1);
		return (
			<div>
				<NavBarNew path={path} />
				<RouteHandler data={data} suggestions={suggestions}/>
				<Database
					show={path.indexOf('variants') === 0}
					suggestions={suggestions}
					data={data}/>
	            <Footer />
            </div>
		);
	}
});

var routes = (
	<Route handler={Application}>
		<DefaultRoute handler={Home}/>
		<Route path='about/:page' handler={About}/>
		<Route path='help' handler={Help}/>
		<Route path='variants' />
		<Route path='variant/:id' handler={VariantDetail}/>
	</Route>
);

var main = document.getElementById('main');

run(routes, HistoryLocation, (Root) => {
  React.render(<Root/>, main);
});
