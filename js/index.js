/*eslint-env browser */
/*global require: false */
'use strict';

var React = require('react');
require('bootstrap/dist/css/bootstrap.css');
var Rx = require('rx');
var vcf = require('vcf.js');
require('rx.binding');
require('rx-dom');
require('../custom.css');

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
	componentDidMount: function () {
		Rx.DOM.click(loadbtn).subscribe(() => {
			this.setState({dialog: true});
		});
	},
	hideDialog: function () {
		this.setState({dialog: false});
	},
	loadData: function (ev) {
		// XXX add error handling
		this.setState({data: vcf.parser()(ev.currentTarget.result)});
	},
	getInitialState: function () {
		return {
			dialog: false,
			data: undefined
		};
	},
	render: function () {
		var {dialog, data} = this.state;
		return (
			<div style={{position: "relative", height: "100px"}}>
				{dialog ? <VCFUpload dataReady={this.loadData} onRequestHide={this.hideDialog}/> : ''}
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
        	<div class="text-center">
        		<h1>BRCA Challenge</h1>
        	</div>
        	)
    }
});

var NavBar = React.createClass({
    render: function() {
    	var {onAbout, onHome, onDB, onMV} = this.props;
        return (
            <ButtonGroup>
                <Button onClick={onHome}>Home</Button>
                <Button onClick={onAbout}>About</Button>
                <Button onClick={onDB}>Database</Button>
                <Button onClick={onMV}>My Variant</Button>
            </ButtonGroup>
        )
    }
});


var Application = React.createClass({
	getInitialState: function () {
		return {about: false, home: false, database: false, myVariant: false};
	},

	onHome: function () {
		this.setState({home: true, about:false, database: false, myVariant: false});
	},

	onAbout: function () {
		this.setState({about: true, home:false, database: false, myVariant: false});
	},

	onDB: function () {
		this.setState({home: false, about:false, database: true, myVariant: false});
	},

	onMV: function () {
		this.setState({home: false, about:false, database: false, myVariant: true});
	},


	render: function () {
		var {about, home, database, myVariant} = this.state;
		return (
			<div>
				<Title />
				<NavBar onAbout={this.onAbout} onHome={this.onHome} 
				onDB={this.onDB} onMV={this.onMV} />
				{about ? <div>Hello</div> : ''}
				{home ? <div><input>variant search</input></div> : ''}
				{database ? <div>Hello db</div> : ''}
				{myVariant ? <div>Hello my variant</div> : ''}				
				<TableView />
			</div>
		);
	}
});

var main = document.getElementById('main');
React.render(<Application />, main);
