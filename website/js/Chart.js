'use strict';

var React = require('react'),
    Highcharts = require('highcharts'),
    _ = require('lodash');

var defaultOptions = {
    credits: {
        enabled: false
    }
};

class Chart extends React.Component {
    // When the DOM is ready, create the chart.
    componentDidMount() {
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
    }
    //Destroy chart before unmount.
    componentWillUnmount() {
        if (this.chart && this.chart.destroy) this.chart.destroy();
    }
    //Create the div which the chart will be rendered to.
    render() {
        return <div id={this.props.container} />;
    }
}

export default Chart;
