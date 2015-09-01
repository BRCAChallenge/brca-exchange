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
var brca_logo = require('./img/brca_logo.png')
var white_bg = require("./img/1px_white.jpg")

var databaseUrl = require('file!../../brca-database.vcf');

var {Col, Row, Input, Modal, Button, ButtonGroup, Navbar, CollapsableNav, Nav, NavItem, 
	Carousel, CarouselItem} = require('react-bootstrap');

var VariantTable = require('./VariantTable');

var TableView = React.createClass({
	render: function () {
		var {data} = this.props;
		return (
			<div style={{position: "relative", height: "100px"}}>
				{data?
					<Row>
						<Col md={10} mdOffset={1}>
							<VariantTable data={data}/>
						</Col>
					</Row>
					: ''}
			</div>
		);
	}
});

var NavBarNew = React.createClass({
    render: function () {
        var {activeButton} = this.props;
        return (
            <div>
            	{/*<img style={{height: 40, width: 40, float: "left"}} src={brca_logo}></img>
                <Navbar brand={<a style={{fontSize: 25, color: "#FF3399"}} href="http://brcaexchange.cloudapp.net">
                			    BRCA Exchange</a>} toggleNavKey={0}> */}
                <Navbar brand={<a href="http://brcaexchange.cloudapp.net">
                			   <img style={{height: 30, width: 30, float: "left", marginBottom: 10}} src={brca_logo} />
                			   <span style={{fontSize: 25, color: "#FF3399"}}>&nbsp;BRCA Exchange&nbsp;&nbsp;&nbsp;</span></a>}
                			   toggleNavKey={0}>

                	<CollapsableNav eventKey={0}>
	                    <Nav navbar>
	                        <NavItem eventKey={1} onClick={() => activeButton('home')}>Home</NavItem>
	                        <NavItem eventKey={2} onClick={() => activeButton('about')}>About</NavItem>
	                        <NavItem eventKey={3} onClick={() => activeButton('database')}>Database</NavItem>
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
            </div>
        )
    }
});

var Title = React.createClass({
    render: function() {
        return (
			<Row style={{marginTop: 30, marginBottom: 30}}>
				<div className="text-center" style={{fontSize: 20, color: "#FF3399"}}>
					<h1>BRCA Challenge</h1>
				</div>
			</Row>
        )
    }
});

var NavBar = React.createClass({
    render: function() {
    	var {activeButton} = this.props;
        return (
        	<Row className="text-center">
	            <ButtonGroup>
	                <Button className="btn-custom" 
	                		onClick={() => activeButton('home')}>Home</Button>
	                <Button className="btn-custom" 
	                		onClick={() => activeButton('about')}>About</Button>
	                <Button className="btn-custom" 
	                		onClick={() => activeButton('database')}>Database</Button>
	                <Button className="btn-custom" 
	                		onClick={() => activeButton('myVariant')}>My Variant</Button>
	            </ButtonGroup>
            </Row>
        )
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
		return(
			<div>
				<div>
					<Row style={{marginTop: 100}}>
						<div className="text-center">place holder for home</div>
					</Row>
				</div>
			</div>
		)
	}
});

var About = React.createClass({
	render: function() {
		return(
			<div>
				<Row style={{marginTop: 100}}>
					<div className="text-center">
						<span>place holder for about</span>
					</div>
				</Row>
			</div>
		)
	}
});

// sketch of function to filter rows on exact matches
function filterData(data, str) {
	var {records, header} = data;
	var filteredRecords = _.filter(data, row => {
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

var MyVariant = React.createClass({
	getInitialState: function () {
		return {
			data: null
		};
	},

	render: function() {
		var {data} = this.state;
		var {show} = this.props;
		return(
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
		)
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
})

var Application = React.createClass({
	getInitialState: function () {
		return {data: null, buttonName: 'null'};
	},
	activeButton: function (buttonName) {
		this.setState({buttonName: buttonName})
	},
	
	hideDialog: function () {
		this.setState({buttonName: "null"});
	},
	render: function () {
		var {buttonName, data} = this.state;
		return (
			<div>
                <NavBarNew activeButton={this.activeButton}/>
				{buttonName === 'about' ? <About /> : ''}
				{buttonName === 'home' ? <Home /> : ''}
				<Database show={buttonName === 'database'} data={data}/>
			</div>
		);
	}
});

var main = document.getElementById('main');
React.render(<Application />, main);
