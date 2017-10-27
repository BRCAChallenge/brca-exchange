'use strict';

var React = require('react'),
    Highcharts = require('highcharts'),
    _ = require('lodash');

var defaultOptions = {
    credits: {
        enabled: false
    }
};

var Chart = React.createClass({
    // When the DOM is ready, create the chart.
    componentDidMount: function () {
        var options = _.merge({}, defaultOptions, this.props.options);
        // Extend Highcharts with modules
        if (this.props.modules) {
            this.props.modules.forEach(function (module) {
                module(Highcharts);
            });
        }
        // Set container which the chart should render to.
        this.chart = new Highcharts[this.props.type || "Chart"](
            this.props.container,
            options
        );
    },
    //Destroy chart before unmount.
    componentWillUnmount: function () {
        this.chart.destroy();
    },
    //Create the div which the chart will be rendered to.
    render: function () {
        return React.createElement('div', { id: this.props.container });
    }
});

module.exports = Chart;
