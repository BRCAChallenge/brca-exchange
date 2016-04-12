/*global module: false, require: false */
'use strict';

var React = require('react');
var _ = require('underscore');
require('muts-needle-plot/src/js/d3-svg-legend');
require('./css/d3Lollipop.css');
var Mutneedles = require("muts-needle-plot");
var PureRenderMixin = require('./PureRenderMixin');

var {Grid, Row, Nav, NavItem} = require('react-bootstrap');


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
        xAxisLabel = 'Coordinate Selection (GRCh38 chr 17)';
        minPos = 43030000;
        maxPos = 43130000;
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
    mixins: [PureRenderMixin],
    render: function () {
        return (
            <div id='brcaLollipop' ref='d3svgBrca'/>
        );
    },
    filterAttributes: function (obj) {
        var oldObj = _(obj).pick('Genomic_Coordinate_hg38', 'Pathogenicity_default');
        var parts = oldObj.Genomic_Coordinate_hg38.split(':');
        var chrCoordinate = parseInt(oldObj.Genomic_Coordinate_hg38.split(':')[parts.length - 2]);
        var alleleChange = _.last(parts);
        var refAllele = alleleChange.split('>')[0];
        var altAllele = alleleChange.split('>')[1];
        if (altAllele.length > refAllele.length) {
            chrCoordinate = String(chrCoordinate) + '-' + String(chrCoordinate + altAllele.length - 1);
        } else {
            chrCoordinate = String(chrCoordinate);
        }
        if (oldObj.Pathogenicity_default == 'Not Yet Classified'){
            oldObj.Pathogenicity_default = "Uncertain";
        }
        if (oldObj.Pathogenicity_default == 'Benign / Little Clinical Significance'){
            oldObj.Pathogenicity_default = "Benign";
        }
        var newObj = {category: oldObj.Pathogenicity_default, coord: chrCoordinate, value: 1, oldData: obj};
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
    getInitialState: function () {
        return {
            brcakey: "BRCA1",
            data: []
        };
    },
    shouldComponentUpdate: function (nextProps, nextState) {
       return (
            //Update if it's the first time rendering
            this.state.data.length === 0 ||
            !_.isEqual(this.props.opts, nextProps.opts) ||
            !_.isEqual(this.state, nextState)
        );
    },
    componentWillReceiveProps: function (newProps) {
        this.fetchData(newProps.opts);
    },
    componentWillMount: function () {
        this.fetchData = _.debounce(this.fetchData, 600, true);
        this.fetchData(this.props.opts);
    },
    fetchData: function (opts) {
        this.props.fetch(opts).subscribe(
            function (d) {
                this.setState({data: d.data});
            }.bind(this));
    },
    onSelect: function (key) {
        this.setState({brcakey: key});
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
