/*global module: false, require: false */
'use strict';

var React = require('react');
var _ = require('underscore');
var {Row, Col, DropdownButton, MenuItem, Grid} = require('react-bootstrap');
require('muts-needle-plot/src/js/d3-svg-legend');
require('./css/d3Lollipop.css');
var Mutneedles = require("muts-needle-plot");
var PureRenderMixin = require('./PureRenderMixin');

var brca12JSON = {
    BRCA1: {
        brcaDomainFile: require('raw!../content/brca1LollipopDomain.json')
    },
    BRCA2: {
        brcaDomainFile: require('raw!../content/brca2LollipopDomain.json')
    }
};
var d3Lollipop = {};

d3Lollipop.drawStuffWithD3 = function(ref, muts, domain, id, varlink, data) {
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
    var config = {variantDetailLink: varlink, minCoord: minPos, maxCoord: maxPos, mutationData: muts, regionData: domain, targetElement: ref.id, legends: legends, colorMap: colorMap };
    var instance =  new Mutneedles(config);
    return function() {
        instance.tip.destroy();
        instance.selectionTip.destroy();
    };
};

var D3Lollipop = React.createClass({
    mixins: [PureRenderMixin],
    render: function () {
        return (
            <div id='brcaLollipop' ref='d3svgBrca'/>
        );
    },
    filterAttributes: function (obj) {
        var oldObj = _(obj).pick('Genomic_Coordinate', 'Clinical_significance');

        var chromosome = oldObj.Genomic_Coordinate.split(':')[1];
        var chrCoordinate = parseInt(oldObj.Genomic_Coordinate.split(':')[1]);
        var alleleChange = oldObj.Genomic_Coordinate.split(':')[2];
        var refAllele = alleleChange.split('>')[0];
        var altAllele = alleleChange.split('>')[1];
        if (altAllele.length > refAllele.length) {
            chrCoordinate = String(chrCoordinate) + '-' + String(chrCoordinate + altAllele.length - 1);
        } else {
            chrCoordinate = String(chrCoordinate);
        }
        if (oldObj.Clinical_significance == '-'){
            oldObj.Clinical_significance = "Unknown";
        }
        var newObj = {category: oldObj.Clinical_significance, coord: chrCoordinate, value: 1, oldData: obj};
        return newObj;
    },
    componentDidMount: function() {
        var {data, brcakey, onRowClick, ...opts} = this.props;
        var subSetData = data.map(this.filterAttributes);
        var d3svgBrcaRef = React.findDOMNode(this.refs.d3svgBrca);
        var domainBRCA = JSON.parse(brca12JSON[brcakey].brcaDomainFile);
        this.cleanupBRCA = d3Lollipop.drawStuffWithD3(d3svgBrcaRef, subSetData, domainBRCA, brcakey, onRowClick);
    },
    componentWillReceiveProps: function(newProps) {
        this.cleanupBRCA();
        var {data, brcakey, onRowClick, ...opts} = newProps;
        var d3svgBrcaRef = React.findDOMNode(this.refs.d3svgBrca);
        while (d3svgBrcaRef.hasChildNodes() ) {
            d3svgBrcaRef.removeChild(d3svgBrcaRef.lastChild);
        }
        var subSetData = data.map(this.filterAttributes);
        var d3svgBrcaRef = React.findDOMNode(this.refs.d3svgBrca);
        while (d3svgBrcaRef.hasChildNodes() ) {
            d3svgBrcaRef.removeChild(d3svgBrcaRef.lastChild);
        }
        var domainBRCA = JSON.parse(brca12JSON[brcakey].brcaDomainFile);
        this.cleanupBRCA = d3Lollipop.drawStuffWithD3(d3svgBrcaRef, subSetData, domainBRCA, brcakey, onRowClick);
    },
    componentWillUnmount: function() {
        this.cleanupBRCA();
    }
});

var Lollipop = React.createClass({
    mixins: [PureRenderMixin],
    getInitialState: function() {
        return {brcakey: "BRCA1"};
    },
    onSelect: function(key) {
        this.setState({brcakey: key});
    },
    render: function () {
        var {data, onHeaderClick, onRowClick, ...opts} = this.props;
        return (
            <Grid>
                <div>
                    <Row style={{marginBottom: '2px', marginTop: '2px'}}>
                        <DropdownButton onSelect={this.onSelect} title="Select Gene" id="bg-vertical-dropdown-1">
                            <MenuItem eventKey="BRCA1">BRCA1</MenuItem>
                            <MenuItem eventKey="BRCA2">BRCA2</MenuItem>
                        </DropdownButton>
                        <span onClick={() => this.props.onHeaderClick('Lollipop Plots')}
                            className='help glyphicon glyphicon-question-sign superscript'/>
                        <D3Lollipop data={this.props.data} key={this.state.brcakey} brcakey={this.state.brcakey} onRowClick={this.props.onRowClick} id='brcaLollipop' ref='d3svgBrca'/>
                    </Row>
                </div>
            </Grid>
        );
    }
});

module.exports = Lollipop;
