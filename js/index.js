/*eslint-env browser */
/*global require: false */
'use strict';

var React = require('react');
require('bootstrap/dist/css/bootstrap.css');
var Rx = require('rx');
var vcf = require('vcf.js');
require('rx.binding');
require('rx-dom');

var Modal = require('react-bootstrap/lib/Modal');
var Input = require('react-bootstrap/lib/Input');
var Row = require('react-bootstrap/lib/Row');
var Col = require('react-bootstrap/lib/Col');

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

var Application = React.createClass({
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

var main = document.getElementById('main');
React.render(<Application />, main);
