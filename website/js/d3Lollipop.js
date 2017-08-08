/*global module: false, require: false */
'use strict';

var React = require('react');
var _ = require('underscore');
require('muts-needle-plot/src/js/d3-svg-legend');
require('./css/d3Lollipop.css');
var Mutneedles = require("muts-needle-plot");
var PureRenderMixin = require('./PureRenderMixin');

var {Grid, Row, Nav, NavItem} = require('react-bootstrap');

//var Spinner = require('spin.js');
//var ReactSpinner = require('react-spinjs');
//var Spinner = require('vitullo-spinner');

var brca12JSON = {
    BRCA1: {
        brcaDomainFile: require('raw!../content/brca1LollipopDomain.json')
    },
    BRCA2: {
        brcaDomainFile: require('raw!../content/brca2LollipopDomain.json')
    }
};


var d3Lollipop = {};

d3Lollipop.drawStuffWithD3 = function(ref, muts, domain, id, varlink) {

    var xAxisLabel = '';
    var minPos = 0;
    var maxPos = 1;
    if (id === 'BRCA1') {
        xAxisLabel = 'Coordinate Selection (GRCh38 chr 17)';
        minPos = 43000000;
        maxPos = 43180000;
    } else if (id === 'BRCA2') {
        xAxisLabel = 'Coordinate Selection (GRCh38 chr 13)';
        minPos = 32300000;
        maxPos = 32410000;
    }
    var legends = {x: xAxisLabel, y: ""};
    var colorMap = {
      // mutation categories
      "Pathogenic": "red",
      "Benign": "lightblue"
    };
    var config = {variantDetailLink: varlink, minCoord: minPos, maxCoord: maxPos, mutationData: muts, regionData: domain, targetElement: ref.id, legends: legends, colorMap: colorMap };
    var instance =  new Mutneedles(config);
    return function() {
        instance.tip.destroy();
        instance.selectionTip.destroy();
    };
};

var D3Lollipop = React.createClass({
    shouldComponentUpdate: () => false,
    render: function () {
        return (
            <div id='brcaLollipop' ref='d3svgBrca'/>
        );
    },
    filterAttributes: function (obj) {
        var oldObj = _(obj).pick('Genomic_Coordinate_hg38', 'Pathogenicity_expert');
        var parts = oldObj.Genomic_Coordinate_hg38.split(':');
        // new format for genomic coordinates, now includes "g.", trim first two characters
        var chrCoordinate = parseInt(parts[1][0] === 'g' ? parts[1].substr(2) : parts[1]);
        var alleleChange = _.last(parts);
        var refAllele = alleleChange.split('>')[0];
        var altAllele = alleleChange.split('>')[1];
        if (altAllele.length > refAllele.length) {
            chrCoordinate = String(chrCoordinate) + '-' + String(chrCoordinate + altAllele.length - 1);
        } else {
            chrCoordinate = String(chrCoordinate);
        }
        if (oldObj["Pathogenicity_expert"] === 'Not Yet Classified') {
            oldObj["Pathogenicity_expert"] = "Uncertain";
        }
        if (oldObj["Pathogenicity_expert"] === 'Benign / Little Clinical Significance') {
            oldObj["Pathogenicity_expert"] = "Benign";
        }
        var newObj = {category: oldObj.Pathogenicity_expert, coord: chrCoordinate, value: 1, oldData: obj};
        return newObj;
    },
    componentWillMount: function() {
        console.log('componentWillmount');
    },
    componentDidMount: function() {
        var {data, brcakey, onRowClick, ...opts} = this.props;
        console.log('componentDidmount start');
        var d3svgBrcaRef = React.findDOMNode(this.refs.d3svgBrca);
        var subSetData = data.map(this.filterAttributes);
        var domainBRCA = JSON.parse(brca12JSON[brcakey].brcaDomainFile);
        // Don't render chart if there's no data recieved yet
        if (this.props.data.length !== 0) {
            this.cleanupBRCA = d3Lollipop.drawStuffWithD3(d3svgBrcaRef, subSetData, domainBRCA, brcakey, onRowClick);
        };
        console.log('componentDidmount end');
    },
    componentWillReceiveProps: function(newProps) {
        // only rebuild plot if number of variants has changed
        console.log('newProps.data.length: ', newProps.data.length);
        console.log('this.props.data.length: ', this.props.data.length);
        console.log('componentWillReveiveProps start');
        if (newProps.data.length !== this.props.data.length) {
            // Don't remove a chart if it wasn't built yet
            if (this.props.data.length !== 0) {
                this.cleanupBRCA();
            };
            var d3svgBrcaRef = React.findDOMNode(this.refs.d3svgBrca);
            var {spinner, spinTarget, data, brcakey, onRowClick, ...opts} = newProps;
            while (d3svgBrcaRef.lastChild) {
                d3svgBrcaRef.removeChild(d3svgBrcaRef.lastChild);
            }
            var subSetData = data.map(this.filterAttributes);
            d3svgBrcaRef = React.findDOMNode(this.refs.d3svgBrca);
            while (d3svgBrcaRef.lastChild ) {
                d3svgBrcaRef.removeChild(d3svgBrcaRef.lastChild);
            }
            var domainBRCA = JSON.parse(brca12JSON[brcakey].brcaDomainFile);
            this.cleanupBRCA = d3Lollipop.drawStuffWithD3(d3svgBrcaRef, subSetData, domainBRCA, brcakey, onRowClick);
        };
        if (newProps.brcakey !== this.props.brcakey) {
            console.log('Inner brcakey changed');
        };
        console.log('componentWillReveiveProps end');
    },
    componentWillUpdate: function () {
        console.log('componentWillUpdate');
    },
    componentDidUpdate: function () {
        console.log('componentDidUpdate');
    },
    componentWillUnmount: function() {
        this.cleanupBRCA();
        console.log('componentWillUnmount');
    }
});

var Lollipop = React.createClass({
    mixins: [PureRenderMixin],
    getInitialState: function () {
        return {
            brcakey: "BRCA1",
            data: [],
        };
    },
    componentWillMount: function () {
        console.log('Outer componentWillMount');
        this.fetchData = _.debounce(this.fetchData, 1000, true);
        this.fetchData(this.props.opts);
    },
    componentDidMount: function () {
        console.log('Outer componentDidMount');
    },
    componentWillReceiveProps: function (newProps) {
        console.log('Outer componentWillReveiveProps');
        this.fetchData(newProps.opts);
    },
    componentWillUpdate: function () {
        console.log('Outer componentWillUpdate');
    },
    componentDidUpdate: function () {
        console.log('Outer componentDidUpdate');
    },
    componentWillUnmount: function() {
        console.log('Outer componentWillUnmount');
    },
    fetchData: function (opts) {
        this.props.fetch(opts).subscribe(
            function (d) {
                this.setState({data: d.data});
                console.log('fetchingData');
            }.bind(this));
    },
    onSelect: function (key) {
        this.setState({brcakey: key});
        console.log('brcakey changed');
    },
    render: function () {
        return (
            <Grid>
                <div>
                    <Row style={{marginBottom: '2px', marginTop: '2px'}}>
                        <Nav bsStyle="tabs" eventKey={0} activeKey={this.state.brcakey} onSelect={this.onSelect} title="Select Gene" id="bg-vertical-dropdown-1">
                            <NavItem eventKey="BRCA1">BRCA1</NavItem>
                            <NavItem eventKey="BRCA2">BRCA2</NavItem>
                        </Nav>
                        <span onClick={() => this.props.onHeaderClick('Lollipop Plots')}/>
                        <D3Lollipop data={this.state.data} opts={this.props.opts} key={this.state.brcakey} brcakey={this.state.brcakey} onRowClick={this.props.onRowClick} id='brcaLollipop' ref='d3svgBrca'/>
                    </Row>
                </div>
            </Grid>
        );
    }
});

module.exports = Lollipop;
