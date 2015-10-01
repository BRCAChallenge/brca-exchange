/*global module: false, require: false */
'use strict';

require('muts-needle-plot/src/js/d3-svg-legend');
var Mutneedles = require("muts-needle-plot");

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
    var legends = {x: xAxisLabel, y: "Variant Count"};
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

module.exports = d3Lollipop;
