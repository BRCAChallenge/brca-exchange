/*global module: false, require: false */
'use strict';

require('muts-needle-plot/src/js/d3-svg-legend');
var mutneedles = require("muts-needle-plot");

var d3Lollipop = {};

if (typeof console === "undefined") {
    window.console = {
        log: function () {}
    };
}

d3Lollipop.drawStuffWithD3 = function(ref, muts, domain) {
    console.log(ref);
    console.log(ref.id);
    console.log('muts: ', muts);
    console.log('domain: ', domain);
    var xAxisLabel = '';
    var minPos = 0;
    var maxPos = 1;
    if (ref.id === 'brca1Lollipop') {
        xAxisLabel = 'BRCA1 Genomic Pos (chr 17)';
        minPos = 41190000;
        maxPos = 41280000;
    } else {
        xAxisLabel = 'BRCA2 Genomic Pos (chr 13)';
        minPos = 32870000;
        maxPos = 32980000;
    }
    var legends = {x: xAxisLabel, y: "Variant Count"};
    var colorMap = {
      // mutation categories
      "Benign": "lightblue",
      "Pathogenic": "red"
    };
    var config = {minCoord: minPos, maxCoord: maxPos, mutationData: muts, regionData: domain, targetElement: ref.id, legends: legends, colorMap: colorMap };
    var instance =  new mutneedles(config);
    return function() {
        instance.tip.destroy();
        instance.selectionTip.destroy();
    };
};

module.exports = d3Lollipop;
