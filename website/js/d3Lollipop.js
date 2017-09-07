/*global module: false, require: false */
'use strict';

var React = require('react');
var _ = require('underscore');
require('muts-needle-plot/src/js/d3-svg-legend');
require('./css/d3Lollipop.css');
var Mutneedles = require("muts-needle-plot");
var PureRenderMixin = require('./PureRenderMixin');

var {Grid, Row, Nav, DropdownButton, MenuItem} = require('react-bootstrap');

var Spinjs = require('spin.js');

var brca12JSON = {
    BRCA1Path: {
        brcaDomainFile: require('raw!../content/brca1LollipopDomain.json')
    },
    BRCA1Allele: {
        brcaDomainFile: require('raw!../content/brca1LollipopDomain.json')
    },
    BRCA2Path: {
        brcaDomainFile: require('raw!../content/brca2LollipopDomain.json')
    },
    BRCA2Allele: {
        brcaDomainFile: require('raw!../content/brca2LollipopDomain.json')
    }
};


var d3Lollipop = {};

d3Lollipop.drawStuffWithD3 = function(ref, muts, domain, brcakey, varlink) {
    var xAxisLabel = '';
    var minPos = 0;
    var maxPos = 1;
    if (brcakey.indexOf('BRCA1') !== -1) {
        xAxisLabel = 'Coordinate Selection (GRCh38 chr 17)';
        minPos = 43000000;
        maxPos = 43180000;
    } else if (brcakey.indexOf('BRCA2') !== -1) {
        xAxisLabel = 'Coordinate Selection (GRCh38 chr 13)';
        minPos = 32300000;
        maxPos = 32410000;
    }
    var legends = {x: xAxisLabel, y: ""};
    var colorMap = {
      // Mutation pathogenicity categories
      "Pathogenic": "red",
      "Benign": "lightblue"
    };
    if (brcakey.indexOf('Allele') !== -1) {
        colorMap = {
            // Allele frequency categories
            "Uncertain": "purple",
            "None": "orange",
            "<0.0001": "red",
            "<0.001": "pink",
            "<0.01": "lightblue",
            "<0.1": "blue",
            "<0.5": "lightgreen",
            ">=0.5": "green"
        };
    };
    var config = {variantDetailLink: varlink, minCoord: minPos, maxCoord: maxPos, mutationData: muts, regionData: domain, targetElement: ref.id, legends: legends, colorMap: colorMap, brcakey: brcakey};
    var instance =  new Mutneedles(config);
    return function() {
        instance.tip.destroy();
        instance.selectionTip.destroy();
    };
};

var D3Lollipop = React.createClass({
    shouldComponentUpdate: () => false,
    getInitialState: function () {
        // loader settings
        var loaderopts = {
          lines: 9, // The number of lines to draw
          length: 9, // The length of each line
          width: 5, // The line thickness
          radius: 14, // The radius of the inner circle
          color: '#EE3124', // #rgb or #rrggbb or array of colors
          speed: 1.9, // Rounds per second
          trail: 40, // Afterglow percentage
          scale: 0.7,
          top: '30%',
          className: 'spinner', // The CSS class to assign to the spinner
        };
        return {
            spinner: new Spinjs(loaderopts),
        };
    },
    render: function () {
        return (
            <div id="spinnerContainer">
                <div id='brcaLollipop' ref='d3svgBrca'/>
            </div>
        );
    },
    startSpinner: function () {
        // Plot target
        var targetElement = document.getElementById('spinnerContainer') || 'spinnerContainer' || document.body;   // Where to append the plot (svg)
        this.state.spinner.spin(targetElement);
    },
    stopSpinner: function () {
        this.state.spinner.stop();
    },
    filterAttributes: function (obj) {
        var oldObj = _(obj).pick('Genomic_Coordinate_hg38', 'Pathogenicity_expert', 'Allele_Frequency');
        var parts = oldObj.Genomic_Coordinate_hg38.split(':');
        // Process genomic coordinates
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
        // Process pathogenicity_expert lableing
        if (oldObj["Pathogenicity_expert"] === 'Not Yet Classified') {
            oldObj["Pathogenicity_expert"] = "Uncertain";
        }
        if (oldObj["Pathogenicity_expert"] === 'Benign / Little Clinical Significance') {
            oldObj["Pathogenicity_expert"] = "Benign";
        }
        // Process Allele Frequency format
        var alleleFreq = parseFloat(oldObj["Allele_Frequency"]);
        var alleleFreqCategory = "";
        if (alleleFreq === 0.0) {
            alleleFreqCategory = "None";
        } else if (alleleFreq < 0.0001) {
            alleleFreqCategory = "<0.0001";
        } else if (alleleFreq < 0.001) {
            alleleFreqCategory = "<0.001";
        } else if (alleleFreq < 0.01) {
            alleleFreqCategory = "<0.01";
        } else if (alleleFreq < 0.1) {
            alleleFreqCategory = "<0.1";
        } else if (alleleFreq < 0.5) {
            alleleFreqCategory = "<0.5";
        } else if (alleleFreq >= 0.5) {
            alleleFreqCategory = ">=0.5";
        } else {
            alleleFreqCategory = "Uncertain";
        };
        var newObj = {category: oldObj.Pathogenicity_expert, coord: chrCoordinate, alleleFreqCategory: alleleFreqCategory, alleleFreq: alleleFreq, value: 1, oldData: obj};
        return newObj;
    },
    componentDidMount: function() {
        this.startSpinner();
        var {dataIsEmpty, data, brcakey, onRowClick, ...opts} = this.props;
        // Redraw the chart when brcakey has changed.
        var d3svgBrcaRef = React.findDOMNode(this.refs.d3svgBrca);
        var subSetData = data.map(this.filterAttributes);
        var domainBRCA = JSON.parse(brca12JSON[brcakey].brcaDomainFile);
        // data.length indicates whether data had PREVIOSLY been loaded.
        // dataIsEmpty indicates if the data is empty AFTER data has been fetched.
        // If 0 data has already been loaded and data will not be empty after data has been fetched,
        //     don't do anything since the chart will be rendered in componentWillReceiveProps.
        //     This occurs when first opening the lollipop chart.
        // If more than 0 data has already been loaded, rerender the chart and stop the spinner.
        //     This occurs when switching gene tabs when there's > 0 matching variants.
        // If the data is still empty after data has been fetched, stop the spinner. Don't render the chart
        //     because there's no data to display.
        //     This occurs when switching gene tabs when there's 0 matching variants.
        if (data.length !== 0) {
            this.cleanupBRCA = d3Lollipop.drawStuffWithD3(d3svgBrcaRef, subSetData, domainBRCA, brcakey, onRowClick);
            // Prevent the spinner from showing up when changing the gene tab
            this.stopSpinner();
        } else if (dataIsEmpty) {
            this.stopSpinner();
        };
    },
    componentWillReceiveProps: function(newProps) {
        // only rebuild plot if number of variants has changed
        if (newProps.data.length !== this.props.data.length) {
            // Don't remove a chart if it wasn't built yet
            if (this.props.data.length !== 0) {
                this.cleanupBRCA();
            };
            var d3svgBrcaRef = React.findDOMNode(this.refs.d3svgBrca);
            var {data, brcakey, onRowClick, ...opts} = newProps;
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
            this.stopSpinner();
        }
        if (newProps.data.length === 0) {
            this.stopSpinner();
        }
    },
    componentWillUnmount: function() {
        // only unmount plot if it was built
        if (typeof this.cleanupBRCA !== "undefined") {
            this.cleanupBRCA();
        };
    }
});

var Lollipop = React.createClass({
    mixins: [PureRenderMixin],
    getInitialState: function () {
        return {
            brcakey: "BRCA1Path",
            data: [],
            dataIsEmpty: false,
        };
    },
    componentWillMount: function () {
        this.fetchData = _.debounce(this.fetchData, 600, true);
        this.fetchData(this.props.opts);
    },
    componentWillReceiveProps: function (newProps) {
        this.fetchData(newProps.opts);
    },
    fetchData: function (opts) {
        this.props.fetch(opts).subscribe(
            function (d) {
                // Added in dataIsEmpty state to make sure that empty data will retrigger componentWillReceiveProps in the child
                //  which is where the stop spinner function is called.
                if (d.data.length === 0) {
                    this.setState({dataIsEmpty: true, data: d.data});
                } else {
                    this.setState({dataIsEmpty: false, data: d.data});
                };
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
                        <Nav bsStyle="tabs" activeKey={this.state.brcakey} onSelect={this.onSelect} id="bg-vertical-dropdown-1">
                            <DropdownButton className={this.state.brcakey.indexOf('BRCA1') > -1 ? 'active' : ''} title="BRCA1" id="brca1-dropdown">
                                <MenuItem eventKey="BRCA1Path">Pathogenicity</MenuItem>
                                <MenuItem eventKey="BRCA1Allele">Allele Frequency</MenuItem>
                            </DropdownButton>
                            <DropdownButton className={this.state.brcakey.indexOf('BRCA2') > -1 ? 'active' : ''} title="BRCA2" id="brca2-dropdown">
                                <MenuItem eventKey="BRCA2Path">Pathogenicity</MenuItem>
                                <MenuItem eventKey="BRCA2Allele">Allele Frequency</MenuItem>
                            </DropdownButton>
                        </Nav>
                        <span onClick={() => this.props.onHeaderClick('Lollipop Plots')}/>
                        <D3Lollipop dataIsEmpty={this.state.dataIsEmpty} data={this.state.data} opts={this.props.opts} key={this.state.brcakey} brcakey={this.state.brcakey} onRowClick={this.props.onRowClick} id='brcaLollipop' ref='d3svgBrca'/>
                    </Row>
                </div>
            </Grid>
        );
    }
});

module.exports = Lollipop;
