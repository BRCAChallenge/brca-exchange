/*global module: false, require: false */
'use strict';

var React = require('react');
var _ = require('underscore');
var {Row, Col, DropdownButton, MenuItem, Grid} = require('react-bootstrap');
require('muts-needle-plot/src/js/d3-svg-legend');
var Mutneedles = require("muts-needle-plot");

var brca12JSON = {
    BRCA1: {
        brcaDomainFile: require('raw!../content/brca1LollipopDomain.json')
    },
    BRCA2: {
        brcaDomainFile: require('raw!../content/brca2LollipopDomain.json')
    }
};
var d3Lollipop = {};

d3Lollipop.drawStuffWithD3 = function(ref, muts, domain, id) {
    var xAxisLabel = '';
    var minPos = 0;
    var maxPos = 1;
    if (id === 'BRCA1') {
        xAxisLabel = 'BRCA1 Genomic Pos (chr 17)';
        minPos = 41190000;
        maxPos = 41280000;
    } else if (id === 'BRCA2') {
        xAxisLabel = 'BRCA2 Genomic Pos (chr 13)';
        minPos = 32880000;
        maxPos = 32980000;
    }
    var legends = {x: xAxisLabel, y: ""};
    var colorMap = {
      // mutation categories
      "Benign": "lightblue",
      "Pathogenic": "red"
    };
    var config = {minCoord: minPos, maxCoord: maxPos, mutationData: muts, regionData: domain, targetElement: ref.id, legends: legends, colorMap: colorMap };
    var instance =  new Mutneedles(config);
    return function() {
        instance.tip.destroy();
        instance.selectionTip.destroy();
    };
};

var D3Lollipop = React.createClass({
    render: function () {
        return (
            <div id='brcaLollipop' ref='d3svgBrca'/>
        );
    },
    filterData: function (obj) {
        if (obj.Gene_symbol === this.props.brcakey && 'Genomic_Coordinate' in obj && 'Clinical_significance' in obj && (obj.Clinical_significance === 'Benign' || obj.Clinical_significance === 'Pathogenic')) {
            return true;
        } else {
            return false;
        }
    },
    filterAttributes: function (obj) {
        var oldObj = _(obj).pick('Genomic_Coordinate', 'Clinical_significance');

        var chrCoordinate = parseInt(oldObj.Genomic_Coordinate.split(':')[1]);
        var refAllele = oldObj.Genomic_Coordinate.split(':')[2].split('>')[0];
        var altAllele = oldObj.Genomic_Coordinate.split(':')[2].split('>')[1];
        if (altAllele.length > refAllele.length) {
            chrCoordinate = String(chrCoordinate) + '-' + String(chrCoordinate + altAllele.length - 1);
        } else {
            chrCoordinate = String(chrCoordinate);
        }
        var newObj = {category: oldObj.Clinical_significance, coord: chrCoordinate, value: 1};
        return newObj;
    },
    componentDidMount: function() {
        var {data, brcakey, ...opts} = this.props;
        var filteredData = data.filter(this.filterData);
        var subSetData = filteredData.map(this.filterAttributes);
        var d3svgBrcaRef = React.findDOMNode(this.refs.d3svgBrca);
        var domainBRCA = JSON.parse(brca12JSON[brcakey].brcaDomainFile);
        this.cleanupBRCA = d3Lollipop.drawStuffWithD3(d3svgBrcaRef, subSetData, domainBRCA, brcakey);
    },
    componentWillRecieveProps: function(newProps) {
        this.setState({data: newProps.data});
        var d3svgBrcaRef = React.findDOMNode(this.refs.d3svgBrca);
        while (d3svgBrcaRef.hasChildNodes() ) {
            d3svgBrcaRef.removeChild(d3svgBrcaRef.lastChild);
        }
    },
    componentWillUpdate: function() {
        var {data, brcakey, ...opts} = this.props;
        var filteredData = data.filter(this.filterData);
        var subSetData = filteredData.map(this.filterAttributes);
        var d3svgBrcaRef = React.findDOMNode(this.refs.d3svgBrca);
        while (d3svgBrcaRef.hasChildNodes() ) {
            d3svgBrcaRef.removeChild(d3svgBrcaRef.lastChild);
        }
        var domainBRCA = JSON.parse(brca12JSON[brcakey].brcaDomainFile);
        this.cleanupBRCA = d3Lollipop.drawStuffWithD3(d3svgBrcaRef, subSetData, domainBRCA, brcakey);
    },
    componentWillUnmount: function() {
        this.cleanupBRCA();
    },
    shouldComponentUpdate: () => true
});

var Lollipop = React.createClass({
    getInitialState: function() {
        return {brcakey: "BRCA1", data: this.props.data};
    },
    onSelect: function(key) {
        this.setState({brcakey: key});
    },
    componentWillReceiveProps: function(newProps) {
        this.setState({data: newProps.data});
    },
    shouldComponentUpdate: () => true,
    render: function () {
        var {data, onHeaderClick, ...opts} = this.props;
        return (
            <Grid>
                <Row>
                    <Col md={8} mdOffset={4}>
                        <h1 id="brca-dna-variant-lollipop">{this.state.brcakey} Lollipop Chart</h1>
                    </Col>
                </Row>
                <div>
                    <DropdownButton onSelect={this.onSelect} title="Select Gene" id="bg-vertical-dropdown-1">
                        <MenuItem eventKey="BRCA1">BRCA1</MenuItem>
                        <MenuItem eventKey="BRCA2">BRCA2</MenuItem>
                    </DropdownButton>
                    <span onClick={() => onHeaderClick('Lollipop Plots')}
                        className='help glyphicon glyphicon-question-sign superscript'/>
                    <D3Lollipop data={this.props.data} key={this.state.brcakey} brcakey={this.state.brcakey} id='brcaLollipop' ref='d3svgBrca'/>
                </div>
            </Grid>
        );
    }
});

module.exports = Lollipop;
