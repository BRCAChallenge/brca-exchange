/*global module: false, require: false */
'use strict';

require('muts-needle-plot/src/js/d3-svg-legend');
var mutneedles = require("muts-needle-plot");
var d3 = require("d3");

var d3Lollipop = {};

d3Lollipop.drawStuffWithD3 = function(ref, data, domain) {
    if (typeof console === "undefined") {
        window.console = {
            log: function () {}
        };
    }
    //var regions = domain;
    //console.log(domain);
    var regions = [
    {"name": "exon1", "coord": "41196311-41197819"},
    {"name": "exon2", "coord": "41199659-41199720"},
    {"name": "exon3", "coord": "41201137-41201211"},
    {"name": "exon4", "coord": "41203079-41203134"},
    {"name": "exon5", "coord": "41209068-41209152"},
    {"name": "exon6", "coord": "41215349-41215390"},
    {"name": "exon7", "coord": "41215890-41215968"},
    {"name": "exon8", "coord": "41219624-41219712"},
    {"name": "exon9", "coord": "41222944-41223255"},
    {"name": "exon10", "coord": "41226347-41226538"},
    {"name": "exon11", "coord": "41228504-41228631"},
    {"name": "exon12", "coord": "41234420-41234592"},
    {"name": "exon13", "coord": "41242960-41243049"},
    {"name": "exon14", "coord": "41243451-41246877"},
    {"name": "exon15", "coord": "41247862-41247939"},
    {"name": "exon16", "coord": "41249260-41249306"},
    {"name": "exon17", "coord": "41251791-41251897"},
    {"name": "exon18", "coord": "41256138-41256278"},
    {"name": "exon19", "coord": "41256884-41256973"},
    {"name": "exon20", "coord": "41258472-41258550"},
    {"name": "exon21", "coord": "41267742-41267796"},
    {"name": "exon22", "coord": "41276033-41276132"},
    {"name": "exon23", "coord": "41277287-41277500"}
];
    var legends = {x: "KRAS genomic pos", y: "Mutation Count"};
    var colorMap = {
      // mutation categories
      "Benign": "lightblue",
      "Pathogenic": "red"
    };

    d3.json('/content/brca1LollipopMuts.json', function(error, data){
        console.log(regions);
        console.log(data);
        var config = {minCoord: 41190000, maxCoord: 41280000, mutationData: data, regionData: regions, targetElement: ref, legends: legends, colorMap: colorMap };
        var instance =  new mutneedles(config);
        return function() {
            instance.tip.destroy();
            instance.selectionTip.destroy();
        };
    });
    //var config = {minCoord: 41190000, maxCoord: 41280000, mutationData: muts, regionData: regions, targetElement: ref, legends: legends, colorMap: colorMap };
    //var config = {minCoord: 25357723, maxCoord: 25403870, mutationData: muts, regionData: regions, targetElement: ref, legends: legends, colorMap: colorMap };
    //var instance =  new mutneedles(config);
    //return function() {
    //    instance.tip.destroy();
    //    instance.selectionTip.destroy();
    //};
};

module.exports = d3Lollipop;
