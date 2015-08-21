/*eslint-env browser */
/*global require: false */
'use strict';

var React = require('react');
require('bootstrap/dist/css/bootstrap.css');
var Rx = require('rx');
var vcf = require('vcf.js');
require('rx.binding');
require('rx-dom');
require('custom.css');

var databaseUrl = require('file!../data/brca.clinvar.vcf');

var {Col, Row, Input, Modal, Button, ButtonGroup} = require('react-bootstrap');

var VariantTable = require('./VariantTable');

// when load is selected, render file upload in pop-up
var loadbtn = document.getElementById('load-vcf');

var VCFUpload = React.createClass({
	getInitialState: function () {
		return {
			file: undefined
		};
	},
	fileChange: function () {
		var {dataReady, onRequestHide} = this.props,
			file = this.refs.file.getInputDOMNode().files[0],
			reader = new FileReader();
		onRequestHide();
		// XXX This timeout allows the UI to update (close dialog) before loading
		// a potentially large file, which will block the UI.
		// This might also be solved by elminating the animation on Modal close,
		// which is probably the source of the problem.
		window.setTimeout(() => {
			reader.onload = dataReady;
			reader.readAsText(file);
		}, 100);
	},
	render: function () {
		var {onRequestHide} = this.props;
		return (
			<Modal title="Import VCF File" onRequestHide={onRequestHide} closeButton={true}>
				<Row>
					<Col md={4} mdOffset={2}>
						<Input ref='file' type='file' onChange={this.fileChange}/>
					</Col>
				</Row>
			</Modal>
		);
	}
});


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
	render: function() {
		return(
			<div>
				<Row style={{marginTop: 100}}>
					<div className="text-center">place holder for home</div>
				</Row>
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
		return (
			<div>
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
		return(
			<div>
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
		var {buttonName,data} = this.state;
		return (
			<div>
				<Title />
				<NavBar activeButton={this.activeButton} />
				{buttonName === 'about' ? <About /> : ''}
				{buttonName === 'home' ? <Home /> : ''}
				{buttonName === 'database' ? <Database data={data}/> : ''}
				{buttonName === 'myVariant' ? <MyVariant dataReady={this.loadData}/> : ''}
			</div>
		);
	}
});

var main = document.getElementById('main');
React.render(<Application />, main);
