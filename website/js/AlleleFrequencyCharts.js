'use strict';

var React = require('react'),
    BarChart = require('./BarChart'),
    util = require('./util'),
    _ = require('lodash');

var alleleFrequencyCharts = function(variant, prop) {
    return <AlleleFrequencyCharts variant={variant} prop={prop} />;
};

class AlleleFrequencyCharts extends React.Component {
    constructor(props) {
        super(props);
        this.state = {
            scaleIndex: 0
        };
    }
    shouldComponentUpdate(nextProps, nextState) {
        let {props, state, refs} = this;
        if (props.variant !== nextProps.variant) {
            return true;
        }

        if (state.scaleIndex !== nextState.scaleIndex) {
            refs[`${props.prop}_alleleFreq2`].getChart().update({
                yAxis: { max: this.scales[nextState.scaleIndex] }
            });
        }

        return false;
    }
    render() {
        let {variant, prop} = this.props;
        let frequencyProps, title, pointFormat;
        if (prop === 'Allele_Frequency_Charts_1000_Genomes') {
            if (!variant['Variant_in_1000_Genomes']) { // eslint-disable-line dot-notation
                return false;
            }
            frequencyProps = [
                {label: 'AFR', prop: 'AFR_Allele_frequency_1000_Genomes'},
                {label: 'AMR', prop: 'AMR_Allele_frequency_1000_Genomes'},
                {label: 'EAS', prop: 'EAS_Allele_frequency_1000_Genomes'},
                {label: 'EUR', prop: 'EUR_Allele_frequency_1000_Genomes'},
                {label: 'SAS', prop: 'SAS_Allele_frequency_1000_Genomes'},
            ];
            title = '1000 Genomes';
            pointFormat =  "{point.y}";
        }
        else if (prop === 'Allele_Frequency_Charts_ExAC') {
            if (!variant['Variant_in_ExAC']) {
                return false;
            }
            frequencyProps = [
                {label: 'AFR', prop: 'Allele_frequency_AFR_ExAC'},
                {label: 'AMR', prop: 'Allele_frequency_AMR_ExAC'},
                {label: 'EAS', prop: 'Allele_frequency_EAS_ExAC'},
                {label: 'FIN', prop: 'Allele_frequency_FIN_ExAC'},
                {label: 'NFE', prop: 'Allele_frequency_NFE_ExAC'},
                {label: 'OTH', prop: 'Allele_frequency_OTH_ExAC'},
                {label: 'SAS', prop: 'Allele_frequency_SAS_ExAC'},
            ];

            title = 'ExAC';
            pointFormat =  "{point.y}<br /><em>({point.count} of {point.number})</em>";
        } else {
            return false;
        }
        // extract all frequencies specified in frequencyProps,
        // if they are all empty fields, do not display chart
        if (_.every(_.values(_.pick(variant, _.map(frequencyProps, e => e.prop))), util.isEmptyField)) {
            return false;
        }

        let categories = [];
        let data = [];

        for (let frequencyType of frequencyProps) {
            if (!util.isEmptyField(variant[frequencyType.prop])) {
                categories.push(frequencyType.label);
                data.push({
                    y: parseFloat(variant[frequencyType.prop]),
                    number: variant[frequencyType.prop.replace("frequency", "number")],
                    count: variant[frequencyType.prop.replace("frequency", "count")]
                });
            }
        }

        // find max frequency of those present, parsing floats from strings and suppressing those that don't parse
        let chart2Max = Math.max.apply(null, _.filter(_.map(_.values(_.pick(variant, _.map(frequencyProps, e => e.prop))), parseFloat), e => isFinite(e)));
        this.scales = [chart2Max, 0.01, 0.001];

        let lineColor = '#D00';
        let fullscaleChartOptions = {
            title: { text: title},
            legend: { enabled: false },
            plotOptions: {
                column: {
                    tooltip: { pointFormat: pointFormat },
                    dataLabels: { enabled: false }
                }
            },
            chart: { spacing: [ 8, 4, 4, 8 ], height: 300 },
            yAxis: {
                max: 1.0,
                tickInterval: null,
                labels: {
                    rotation: 45,
                    style: { fontSize: "8px" },
                    x: -2,
                    y: 2
                },
                plotLines: [{
                    color: lineColor,
                    dashStyle: 'Dash',
                    value: 0.01,
                    width: 1,
                    zIndex: 5,
                    label: {
                        text: '1%',
                        align: 'left',
                        style: { color: lineColor }
                    }
                }]
            },
            xAxis: { categories: categories },
            series: [{ data: data, allowPointSelect: false }]
        };

        let scaledChartOptions = {
            title: { text: `${title} (scaled)`},
            legend: { enabled: false },
            plotOptions: {
                column: {
                    tooltip: { pointFormat: pointFormat },
                    dataLabels: { enabled: false }
                }
            },
            chart: { spacing: [ 8, 10, 4, 0 ], height: 300 },
            yAxis: {
                max: chart2Max,
                tickInterval: null,
                labels: {
                    rotation: 45,
                    style: { fontSize: "8px" },
                    x: -2,
                    y: 2
                },
                plotLines: [{
                    color: lineColor,
                    dashStyle: 'Dash',
                    value: 0.01,
                    width: 1,
                    zIndex: 5,
                    label: {
                        text: '1%',
                        align: 'left',
                        style: { color: lineColor }
                    }
                }]
            },
            xAxis: { categories: categories },
            series: [{ data: data, allowPointSelect: false }]
        };

        return (
            <div>
                <div className='alleleFrequencyChart'>
                    <BarChart ref={`${prop}_alleleFreq2`} container={`${prop}_alleleFreq1`} options={fullscaleChartOptions} onClick={() => alert ('foo')}/>
                </div>
                <div className='alleleFrequencyChart' onClick={this.toggleScale.bind(this, `${prop}_alleleFreq2`)}>
                    <BarChart ref={`${prop}_alleleFreq2`} container={`${prop}_alleleFreq2`} options={scaledChartOptions} />
                    <div style={{textAlign: 'center', color: 'grey', fontSize: '12px'}}>(click chart to change scale)</div>
                </div>
            </div>
        );
    }

    toggleScale() {
        let idx = this.state.scaleIndex;
        this.setState({
            scaleIndex: (idx + 1) % this.scales.length
        });
    }
}

module.exports = alleleFrequencyCharts;
