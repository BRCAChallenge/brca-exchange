'use strict';

var React = require('react'),
    _ = require('lodash'),
    Chart = require('./Chart');

var defaultOptions = {
    chart: {
            plotBackgroundColor: null,
            plotBorderWidth: null,
            plotShadow: false,
            type: 'pie'
    },
    tooltip: {
        pointFormat: '{series.name}: <b>{point.percentage:.1f}%</b>'
    },
    plotOptions: {
        pie: {
        allowPointSelect: true,
        cursor: 'pointer',
            dataLabels: {
                enabled: false
            },
            showInLegend: true
        }
    }
};

var PieChart = React.createClass({
    render() {
        var options = _.merge({}, defaultOptions, this.props.options);
        return <Chart ref={"chart"} container={this.props.container} options={options} />;
    },
    getChart() {
        return this.refs.chart.chart;
    }
});

module.exports = PieChart;
