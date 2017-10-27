'use strict';

var React = require('react'),
    _ = require('lodash'),
    Highcharts = require('highcharts'),
    Chart = require('./Chart');

require('highcharts/modules/broken-axis')(Highcharts);

var defaultOptions = {
    chart: {
        plotBackgroundColor: null,
        plotBorderWidth: null,
        plotShadow: false,
        type: 'column'
    },
    plotOptions: {
        column: {
            allowPointSelect: true,
            cursor: 'pointer',
            dataLabels: {
                enabled: true,
                format: '<b>{point.y}</b>',
                style: {
                    color: (Highcharts.theme && Highcharts.theme.contrastTextColor) || 'black'
                }
            }
        }
    },
    yAxis: {
        title: false,
        tickInterval: 1000,
        events: {
            pointBreak: pointBreakColumn
        }
    }
};

Highcharts.wrap(Highcharts.Axis.prototype, 'getLinePath', function (proceed, lineWidth) {
    var axis = this,
        path = proceed.call(this, lineWidth),
        x = path[1],
        y = path[2];

    Highcharts.each(this.breakArray || [], function (brk) {
        if (axis.horiz) {
            x = axis.toPixels(brk.from);
            path.splice(3, 0,
                'L', x - 4, y, // stop
                'M', x - 9, y + 5, 'L', x + 1, y - 5, // left slanted line
                'M', x - 1, y + 5, 'L', x + 9, y - 5, // higher slanted line
                'M', x + 4, y
            );
        } else {
            y = axis.toPixels(brk.from);
            path.splice(3, 0,
                'L', x, y - 4, // stop
                'M', x + 5, y - 9, 'L', x - 5, y + 1, // lower slanted line
                'M', x + 5, y - 1, 'L', x - 5, y + 9, // higher slanted line
                'M', x, y + 4
            );
        }
    });
    return path;
});

/**
 * On top of each column, draw a zigzag line where the axis break is.
 */
function pointBreakColumn(e) {
    var point = e.point,
        brk = e.brk,
        shapeArgs = point.shapeArgs,
        x = shapeArgs.x,
        y = this.translate(brk.from, 0, 1, 0, 1),
        w = shapeArgs.width,
        key = ['brk', brk.from, brk.to],
        path = ['M', x, y, 'L', x + w * 0.25, y + 4, 'L', x + w * 0.75, y - 4, 'L', x + w, y];

    if (!point[key]) {
        point[key] = this.chart.renderer.path(path)
            .attr({
                'stroke-width': 2,
                stroke: point.series.options.borderColor
            })
            .add(point.graphic.parentGroup);
    } else {
        point[key].attr({
            d: path
        });
    }
}

var BarChart = React.createClass({
    render() {
        var options = _.merge({}, defaultOptions, this.props.options);
        return <Chart ref={"chart"} container={this.props.container} options={options} />;
    },
    getChart() {
        return this.refs.chart.chart;
    }
});

module.exports = BarChart;
