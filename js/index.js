/*eslint-env browser */
/*global require: false */
'use strict';

var React = require('react');
require('bootstrap/dist/css/bootstrap.css');
var Rx = require('rx');
var vcf = require('vcf.js');
require('rx.binding');
require('rx-dom');
require('css/custom.css');
var _ = require('underscore');
var brca_logo = require('./img/brca_logo.png');
var ga4gh_logo = require('./img/ga4gh-logo-less.png');

var Markdown = require('react-remarkable');
var content = {
	home: require('../content/home.md'),
	'about-history': require('../content/history.md'),
	'about-what': require('../content/brca1_2.md'),
	'about-variation': require('../content/variationAndCancer.md')
};


var databaseUrl = require('file!../../brca-database.vcf');

var {Col, Row, Input, Button, Navbar, CollapsableNav, Nav,
	NavItem, DropdownButton, MenuItem} = require('react-bootstrap');


var VariantTable = require('./VariantTable');

var NavBarNew = React.createClass({
	close: function () {
		this.setState({open: false});
	},
	render: function () {
		var {activeButton} = this.props;
		return (
			<Navbar>
				<a className="navbar-brand" href="#">
					<img style={{height: 20, width: 20, display: 'inline-block'}} src={brca_logo} alt="brca logo"/>
					BRCA Exchange
				</a>
				<CollapsableNav>
					<Nav navbar>
						<NavItem onClick={() => activeButton('home')}>Home</NavItem>
						<DropdownButton onSelect={this.close} title='About'>
							<MenuItem onClick={() => activeButton('about-history')}>
								History of the BRCA Exchange
							</MenuItem>
							<MenuItem onClick={() => activeButton('about-what')}>
								What are BRCA1 and BRCA2?
							</MenuItem>
							<MenuItem onClick={() => activeButton('about-variation')}>
								BRCA Variation and Cancer
							</MenuItem>
						</DropdownButton>
						<NavItem onClick={() => activeButton('database')}>Database</NavItem>
					</Nav>
					<Nav navbar right>
						<NavItem href='#'><input placeholder="Search Variant"></input>
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
			<div>
				<Row style={{marginTop: 100}}>
					<div className="text-center">
						<Markdown options={{html: true}} source={content.home} />
					</div>
				</Row>
				<div className="logo-block">
					<div className="logos">
						<ul>
							<li><img src={ga4gh_logo} alt="ga4gh logo" /></li>
						</ul>
					</div>
				</div>

			</div>
		);
	}
});

var About = React.createClass({
	render: function() {
		var {contentKey} = this.props;

		return (
			<div>
				<Row style={{marginTop: 100}}>
					<div className="text-center">
						<Markdown options={{html: true}} source={content[contentKey]} />
					</div>
				</Row>
			</div>
		);
	}
});

// sketch of function to filter rows on exact matches
function filterData(data, str) { //eslint-disable-line no-unused-vars
	var {records, header} = data;
	var filteredRecords = _.filter(records, row => {
		// row = {
		//   chrom: "17",
		//   pos: 1234,
		//   hgvs: "NC_0001:1234T>C"
		// }
		//
		return _.find(_.values(row), col => col.indexOf(str) !== -1);

	});
	return {
		records: filteredRecords,
		header: header
	};
}

var Database = React.createClass({
	getInitialState: function () {
		return {
			data: null
		};
	},
	componentWillMount: function (){
		Rx.DOM.get(databaseUrl).subscribe(data => this.setState({data: vcf.parser()(data.response)}));
	},

	render: function () {
		var {data} = this.state;
		var {show} = this.props;
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
							<input><span className="glyphicon glyphicon-search"></span></input>
						</div>
					</Row>
					<Row className="text-center" style={{fontSize: 12, color: "grey"}}>
						search for known variants
					</Row>
				</div>

				<div style={{position: "relative", height: "100px"}}>
					{data?
						<Row>
							<Col md={10} mdOffset={1}>
								<VariantTable data={data}/>
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
					{data?
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

var startsWith = (pat, str) => str && str.indexOf(pat) === 0;

var Application = React.createClass({
	getInitialState: function () {
		return {data: null, buttonName: null};
	},
	activeButton: function (buttonName) {
		this.setState({buttonName: buttonName});
	},
	render: function () {
		var {buttonName, data} = this.state;
		return (
			<div>
                <NavBarNew activeButton={this.activeButton}/>
				{startsWith('about', buttonName) ? <About contentKey={buttonName} /> : ''}
				{buttonName === 'home' ? <Home /> : ''}
				<Database show={buttonName === 'database'} data={data}/>
			</div>
		);
	}
});

var main = document.getElementById('main');
React.render(<Application />, main);
