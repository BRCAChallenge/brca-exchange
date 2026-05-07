'use strict';

import React from 'react';
import _ from 'lodash';
import Chart from './Chart';

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

class PieChart extends React.Component {
    constructor(props) {
	super(props);
	this.chartRef = React.createRef();
    }
    render() {
        var options = _.merge({}, defaultOptions, this.props.options);
        return <Chart ref={this.chartRef} container={this.props.container} options={options} />;
    }
    getChart() {
	return this.chartRef.current?.chart;
    }
}

export default PieChart;
