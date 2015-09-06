/*eslint-env browser */
/*global require: false */
'use strict';

var React = require('react');
require('bootstrap/dist/css/bootstrap.css');
require('font-awesome-webpack');
var Rx = require('rx');
var vcf = require('vcf.js');
require('rx.binding');
require('rx-dom');
require('css/custom.css');
var _ = require('underscore');
//var brcaLogo = require('./img/brca_logo.png');
var ga4ghLogo = require('./img/ga4gh-logo-less.png');
var brcaLogoWithText = require('./img/BRCA-logo-with-text.png');
var hvpLogo = require('./img/hvp_logo.png');
var UNESCOLogo = require('./img/UNESCO-logo.jpg');
var ENIGMALogo = require('./img/enigma_logo.png');
var CIMBALogo = require('./img/cimba_logo.png');
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
			<div {...otherProps} dangerouslySetInnerHTML={{__html: html}} />
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
	render: function () {
		return (
			<Navbar toggleNavKey={0}>
				<a className="navbar-brand" href="http://brcaexchange.org">
					<span style={{margin: 10, fontSize: 30, color: "#FF3399"}}>
						BRCA Exchange
					</span>
				</a>
				<CollapsableNav eventKey={0}>
					<Nav navbar>
						<NavLink eventKey={1} to='/'>Home</NavLink>
						<DropdownButton eventKey={2} ref='about' title='About'>
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
						<NavLink eventKey={3} to='/variants'>Variants</NavLink>
						<NavLink to='/help'>Help</NavLink>
					</Nav>
					<Nav navbar right>
						<NavItem eventKey={1} href='#'><input placeholder="Search Variant"></input>
							<Button className='btn-xs' style={{border: 0}}>
								<span className="glyphicon glyphicon-search"></span>
							</Button>
						</NavItem>
					</Nav>
				</CollapsableNav>
			</Navbar>
		);
	}
});

var Home = React.createClass({
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

	render: function() {
		return (
			<Grid>
				<Row style={{marginTop: 10}}>
					<Col md={4} mdOffset={4}>
						<Well>
							<div className="text-center">
								<input placeholder="Search Variant">
									<span className="glyphicon glyphicon-search"/>
								</input>
								<span className="glyphicon glyphicon-question-sign superscript"/>
							</div>
						</Well>
					</Col>
				</Row>
				<Row style={{marginTop: 100}}>
					<Col md={8} mdOffset={2}>
						<RawHTML html={content.pages.home} />
					</Col>
				</Row>
				<Row className='logo-block'>
					<Col md={6} mdOffset={3}>
						<ul className='logos'>
							<li><a href="http://genomicsandhealth.org">
								<img src={ga4ghLogo} alt="ga4gh logo" />
							</a></li>
							<li><a href="http://brcaexchange.org">
								<img src={brcaLogoWithText} alt="brca exchange logo" />
							</a></li>
							<li><a href="http://www.humanvariomeproject.org">
								<img src={hvpLogo} alt="human variome project logo" />
							</a></li>
							<li><a href="http://unesco.org">
								<img src={UNESCOLogo} alt="UNESCO logo" />
							</a></li>
							<br></br>
							<li><a href="http://enigmaconsortium.org">
								<img src={ENIGMALogo} alt="ENIGMA logo" />
							</a></li>
							<li><a href="http://apps.ccge.medschl.cam.ac.uk/consortia/cimba//">
								<img src={CIMBALogo} alt="CIMBA logo" />
							</a></li>
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
				<Row style={{marginTop: 100}}>
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
				<Row style={{marginTop: 100}}>
					<Col md={8} mdOffset={2}>
						<RawHTML ref='content' html={content.pages.help} />
					</Col>
				</Row>
			</Grid>
		);
	}
});

var Database = React.createClass({
	mixins: [Navigation],
	showVariant: function (id) {
		this.transitionTo(`/variant/${id}`);
	},
	showHelp: function (title) {
		this.transitionTo(`/help#${slugify(title)}`);
	},
	render: function () {
		var {show, data} = this.props;
		return (
			<div style={{display: show ? 'block' : 'none'}}>
				<div>
					<Row style={{marginTop: 50, marginBottom: 50}}>
						<div className="text-center">
							<span>place holder for database summary</span>
						</div>
					</Row>
					<Row style={{marginTop: 10}}>
						<div className="text-center">
							<input><span className="glyphicon glyphicon-search"/></input>
							<span className="glyphicon glyphicon-question-sign superscript"/>
						</div>
					</Row>
					<Row className="text-center" style={{fontSize: 12, color: "grey"}}>
						search for known variants
					</Row>
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
							<Button>Download</Button>
						</Col>
					</Row>
				</div>


				<div style={{position: "relative", height: "100px"}}>
					{data ?
						<Row>
							<Col md={10} mdOffset={1}>
								<VariantTable
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
	render: function() {
		var {data, params: {id}} = this.props,
			variant = (data && data.records[id]) || {};

		variant = _.omit(variant, ['__HEADER__']);
		var rows = _.map(variant, (v, k) => <tr key={k}><td>{k}</td><td>{v}</td></tr>);
		return (
			<Grid>
				<Row style={{marginTop: 100}}>
					<div className='text-center'>
						<h3>Variant Detail</h3>
						<p>Placeholder for a better variant details page</p>
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
				<Database show={path === 'variants'} data={data}/>
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
