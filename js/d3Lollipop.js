// d3Lollipop.js

require('muts-needle-plot/src/js/d3-svg-legend');
var mutneedles = require("muts-needle-plot");


var d3Lollipop = {};

d3Lollipop.drawStuffWithD3 = function(ref, data) {
    
    var muts = [
      {
        "coord":"25398273-25398282",
        "category": "cluster1",
        "value":461
      },
      {
        "coord":"25378645-25378645",
        "category": "cluster2",
        "value":4
      },
      {
        "coord":"25380257-25380257",
        "category": "cluster3",
        "value":36
      }
    ];
    var regions = [
      {"name": "exon1", "coord": "25403685-25403865"},
      {"name": "exon2", "coord": "25398208-25398329"},
      {"name": "exon3", "coord": "25380168-25380346"},
      {"name": "exon4", "coord": "25378548-25378707"},
      {"name": "exon5", "coord": "25357723-25362845"}
    ];
    var legends = {x: "KRAS genomic pos", y: "Mutation Count"}
    var colorMap = {
      // mutation categories
      "missense_variant": "yellow",
      "frameshift_variant": "blue",
      "stop_gained": "red",
      "stop_lost": "orange",
      "synonymous_variant": "lightblue"
    };


    var config = {minCoord: 25357723, maxCoord: 25403870, mutationData: muts, regionData: regions, targetElement: ref, legends: legends, colorMap: colorMap }

    var instance =  new mutneedles(config);

    return function() {
        instance.tip.destroy();
        instance.selectionTip.destroy();
    };
}

module.exports = d3Lollipop;
