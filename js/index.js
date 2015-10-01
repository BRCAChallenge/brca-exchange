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

var brcaLogo = require('./img/BRCA-Exchange-tall-tranparent.png');
var logos = require('./logos');
var slugify = require('./slugify');

var content = require('./content');
var brca12JSON = {
    BRCA1: {
        brcaMutsFile: require('raw!../content/brca1LollipopMuts.json'),
        brcaDomainFile: require('raw!../content/brca1LollipopDomain.json')
    },
    BRCA2: {
        brcaMutsFile: require('raw!../content/brca2LollipopMuts.json'),
        brcaDomainFile: require('raw!../content/brca2LollipopDomain.json')
    }
};

var databaseUrl = require('../../enigma-database.tsv');
var databaseKey = require('../databaseKey');

var {Grid, Col, Row, Input, Navbar, Nav, Table,
	DropdownButton, MenuItem, Modal, Button} = require('react-bootstrap');


var VariantTable = require('./VariantTable');
var VariantSearch = require('./VariantSearch');
var {Navigation, State, Link, Route, RouteHandler,
	HistoryLocation, run, DefaultRoute} = require('react-router');

var navbarHeight = 70; // XXX This value MUST match the setting in custom.css

var d3Lollipop = require('./d3Lollipop');


var variantPathJoin = row => _.map(databaseKey, k => encodeURIComponent(row[k])).join('@@');
var variantPathSplit = id => _.object(databaseKey, _.map(id.split(/@@/), decodeURIComponent));

if (typeof console === "undefined") {
    window.console = {
        log: function () {}
    };
}

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
						<NavLink onClick={this.close} to='/about/lollipop'>
							DNA Variant BRCA Lollipop Plots
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
            <div className="container footer">
                <div className="col-sm-5 left-footer">
                    <ul>
                        <li><a href="/home">Home</a></li>
                        <li><a href="/about/history">About</a></li>
                        <li><a href="/variants">Variants</a></li>
                        <li><a href="/help">Help</a></li>
                    </ul>
                </div>
                <div className="col-sm-2 logo-footer">
                    <img href="#" src={brcaLogo} alt="brca exchange logo" />
                </div>
                <div className="col-sm-5 right-footer">
                    <ul>
                        <li><DisclaimerModal /></li>
                        <li><a href="mailto:brca-exchange-contact@genomicsandhealth.org?subject=BRCA Exchange website">contact us</a></li>
                        <li>
                            <a href="https://github.com/BD2KGenomics/brca-website">
                                source code
                            </a>
                        </li>
                    </ul>
                </div>
            </div>
        );
    }
});

var DisclaimerModal = React.createClass({
    getInitialState() {
        return { showModal: false };
    },
    close() {
        this.setState({ showModal: false });
    },
    open() {
        this.setState({ showModal: true });
    },
    onRequestHide() {
        this.setState({ showModal: false });
    },
    render() {
        return (
            <div style={{display: "inline"}}>
                <a onClick={this.open}>disclaimer</a>
                {this.state.showModal ?
                    <Modal onHide={this.close}>
                        <RawHTML html={content.pages.disclaimer} />
                        <div className = "close-button">
                            <Button onClick={this.close}>close</Button>
                        </div>
                    </Modal> : null }
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
					<Col smOffset={3} sm={6}>
						<VariantSearch
							id='home-search'
							suggestions={suggestions}
							onSearch={this.onSearch}/>
					</Col>
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

var _toSpace = s => s.replace(/_/g, ' ');

var Key = React.createClass({
	render() {
		var {onClick, tableKey} = this.props,
			words = tableKey.replace(/_/g, ' ').split(' ');
		return (
			 <td className='help-target'>
				{words.slice(0, words.length - 1).join(' ')}
				{' '}
				<span className="text-nowrap">
					{words[words.length - 1]}
					<span role='button' onClick={onClick}
						className='help glyphicon glyphicon-question-sign superscript'/>
				</span>
			 </td>
		);
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
				<Key tableKey={k} onClick={() => this.showHelp(_toSpace(k))} />
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
	'Assertion_method_citation',
	'URL'
];

var flatmap = (coll, fn) => _.flatten(_.map(coll, fn), true);
var minSuggestion = 3; // minimum length of string to use in autocomplete
var rowWords = row => flatmap(_.values(_.omit(row, dontSuggest)),
		v => v.toLowerCase().split(/\s+/));

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

var D3Lollipop = React.createClass({
    render: function () {
        return (
            <div id='brcaLollipop' ref='d3svgBrca'/>
        );
    },
    componentDidMount: function() {
        var d3svgBrcaRef = React.findDOMNode(this.refs.d3svgBrca);
        var mutsBRCA = JSON.parse(brca12JSON[this.props.brcakey].brcaMutsFile);
        var domainBRCA = JSON.parse(brca12JSON[this.props.brcakey].brcaDomainFile);
        this.cleanupBRCA = d3Lollipop.drawStuffWithD3(d3svgBrcaRef, mutsBRCA, domainBRCA, this.props.brcakey);
    },
    componentWillUnmount: function() {
        this.cleanupBRCA();
    },
    shouldComponentUpdate: () => false
});

var Lollipop = React.createClass({
    showHelp: function (title) {
        this.transitionTo(`/help#${slugify(title)}`);
    },
    getInitialState: function() {
        return {brcakey: "BRCA1"};
    },
    onSelect: function(key) {
	    this.setState({brcakey: key});
    },
    render: function () {
        return (
            <Grid>
                <Row>
                    <Col md={8} mdOffset={4}>
                        <h1 id="brca-dna-variant-lollipop">{this.state.brcakey} Lollipop Chart</h1>
                    </Col>
                </Row>
                <div>
                    <span onClick={() => this.showHelp('lollipop-plots')}
                        className='help glyphicon glyphicon-question-sign superscript'/>
                    <DropdownButton onSelect={this.onSelect} title="Select Gene" id="bg-vertical-dropdown-1">
                        <MenuItem eventKey="BRCA1">BRCA1</MenuItem>
                        <MenuItem eventKey="BRCA2">BRCA2</MenuItem>
                    </DropdownButton>
                    <D3Lollipop key={this.state.brcakey} brcakey={this.state.brcakey} id='brcaLollipop' ref='d3svgBrca'/>
                </div>
            </Grid>
        );
    }
});

var routes = (
	<Route handler={Application}>
		<DefaultRoute handler={Home}/>
		<Route path='about/lollipop' handler={Lollipop}/>
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
