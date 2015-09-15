/*eslint-env browser */
/*global require: false */
'use strict';

require('./favicons');
var React = require('react');
var PureRenderMixin = require('./PureRenderMixin'); // deep-equals version of PRM
require('bootstrap/dist/css/bootstrap.css');
require('font-awesome-webpack');
var Rx = require('rx');
var vcf = require('vcf.js');
require('rx.binding');
require('rx-dom');
require('css/custom.css');
var _ = require('underscore');
//var brcaLogo = require('./img/brca_logo.png');
var logos = require('./logos');
var slugify = require('./slugify');

var content = require('./content');

var databaseUrl = require('file!../../enigma-database.tsv');

var {Well, Grid, Col, Row, Input, Button, Navbar, CollapsableNav, Nav, Table,
	NavItem, DropdownButton, MenuItem, Panel} = require('react-bootstrap');


var VariantTable = require('./VariantTable');
var {Navigation, State, Link, Route, RouteHandler,
	HistoryLocation, run, DefaultRoute} = require('react-router');

var merge = (...objs) => _.extend({}, ...objs);

// add unique id to variant table
function addId(data) {
	return _.map(data, (r, i) => merge({id: i}, r));
}

function cutTrailingNewLine(string) {
    if (string[string.length - 1] === "\n") {
        return string.slice(0, string.length - 1);
    }
    return string;
}

function readTsv(response) {
	var [header, ...records] = cutTrailingNewLine(response).split("\n");
	var keys = header.split("\t");
    var rows = _.map(records, row => row.split("\t"));
    return {
        records: addId(_.map(rows, row => _.object(keys, row)))
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

var VariantSearch = React.createClass({
	mixins: [Navigation],
	onClick: function () {
		var value = React.findDOMNode(this.refs.input).value;
		this.transitionTo(`/variants?search=${value}`);
	},
	onKeyDown: function (ev) {
		if (ev.key === 'Enter') {
			this.transitionTo(`/variants?search=${ev.target.value}`);
		}
	},
	render: function () {
		return (
			<div className="search-box">
				<input ref='input' onKeyDown={this.onKeyDown} placeholder="Search Variant"></input>
				<Button onClick={this.onClick} className='btn-xs'>
					<span>&nbsp;&nbsp;</span>
                    <span className="glyphicon glyphicon-search"></span>
                    <span onClick={() => this.showHelp('Searching')}
						className="glyphicon glyphicon-question-sign superscript help"/>
				</Button>
			</div>
		);
	}
});

var NavBarNew = React.createClass({
	close: function () {
		this.refs.about.setState({open: false});
	},
	render: function () {
		return (
			<div className="navbar-container">
            <Navbar className="navbar-fixed-top">
				<a className="navbar-brand" href="http://brcaexchange.org">
					<span>
						<b className="BRCA">BRCA</b>
                        <span className="exchange"> Exchange</span> 
					</span>
				</a>
					<Nav navbar right>
						<NavLink to='/'>Home</NavLink>
						<DropdownButton ref='about' title='About'>
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
	showHelp: function (title) {
		this.transitionTo(`/help#${slugify(title)}`);
	},

	render: function() {
		var logoItems = _.map(logos, ({id, logo, url}) => (
			<li key={id}><a href={url}>
				<img src={logo} alt={id + ' logo'} />
			</a></li>
		));
		return (
			<Grid>
				<Row>
				   	<div className='text-center'>
						<VariantSearch />
					</div>
				</Row>
				<Row>
					<Col md={8} mdOffset={2}>
						<RawHTML html={content.pages.home} />
					</Col>
				</Row>
				<Row className='logo-block'>
					<Col md={6} mdOffset={3}>
						<ul className='logos'>
							{logoItems}
						</ul>
					</Col>
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
					window.scrollTo(0, el.getBoundingClientRect().top);
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
	showVariant: function (id) {
		this.transitionTo(`/variant/${id}`);
	},
	showHelp: function (title) {
		this.transitionTo(`/help#${slugify(title)}`);
	},
	createDownload: function (ev) {
		// XXX This is a bit horrible. In order to build the tsv lazily (on button click, instead
		// of on every tabe update), we catch the mousedown event and modify the href on the
		// anchor element, behind the back of react. I don't believe this will cause any
		// problems, but it's something to be aware of if react starts doing something strange.
		// Also needs to be tested cross-browser. We should not offer download on browsers that don't allow
		// client-driven download.
		var data = this.refs.table.getData(),
			keys = _.keys(data[0]),
			tsvRows = _.map(data, obj => _.map(keys, k => obj[k]).join('\t')).join('\n'), // use os-specific line endings?
			tsv = keys.join('\t') + '\n' + tsvRows;
		ev.target.href = URL.createObjectURL(new Blob([tsv], { type: 'text/tsv' }));
	},
	render: function () {
		var {show, data} = this.props,
			{search} = this.getQuery();
		return (
			<div style={{display: show ? 'block' : 'none'}}>
				<div>
					<Row>
						<Col md={4} mdOffset={4}>
							<Panel header='Advanced filtering'>
								<DropdownButton title="Gene">
									<MenuItem>BRCA1</MenuItem>
									<MenuItem>BRCA2</MenuItem>
								</DropdownButton>
								<span className="glyphicon glyphicon-question-sig superscript"/>
								<DropdownButton title="Exon">
									<MenuItem>Any</MenuItem>
									<MenuItem>1</MenuItem>
									<MenuItem>2</MenuItem>
									<MenuItem>3</MenuItem>
									<MenuItem>4</MenuItem>
									<MenuItem>5</MenuItem>
								</DropdownButton>
								<span className="glyphicon glyphicon-question-sign superscript"/>
							</Panel>
						</Col>
						<Col md={1}>
							<Button download="variants.tsv" href="#" onMouseDown={this.createDownload}>Download</Button>
						</Col>
					</Row>
				</div>


				<div style={{position: "relative", height: "100px"}}>
					{data ?
						<Row>
							<Col md={10} mdOffset={1}>
								<VariantTable
									ref='table'
									filterValues={{globalSearch: search || ''}}
									data={data.records}
									onHeaderClick={this.showHelp}
									onRowClick={this.showVariant}/>
							</Col>
						</Row>
						: ''}
				</div>
			</div>
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
			variant = (data && data.records[id]) || {};

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

var Application = React.createClass({
	mixins: [State],
	getInitialState: function () {
		return {data: null};
	},
	componentWillMount: function (){
		Rx.DOM.get(databaseUrl).subscribe(data =>
			this.setState({data: readTsv(data.response)}));
	},
	render: function () {
		var {data} = this.state;
		var path = this.getPath().slice(1);
		return (
			<div>
				<NavBarNew />
				<RouteHandler data={data}/>
				<Database show={path.indexOf('variants') === 0} data={data}/>
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
