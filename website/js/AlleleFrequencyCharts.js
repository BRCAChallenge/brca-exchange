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
        if (prop === 'Allele_Frequency_Charts_Exome_GnomAD') {
            if (!variant['Variant_in_GnomAD']) {
                return false;
            }
            frequencyProps = [
                {label: 'AFR', prop: 'Allele_frequency_exome_AFR_GnomAD'},
                {label: 'AMR', prop: 'Allele_frequency_exome_AMR_GnomAD'},
                {label: 'ASJ', prop: 'Allele_frequency_exome_ASJ_GnomAD'},
                {label: 'EAS', prop: 'Allele_frequency_exome_EAS_GnomAD'},
                {label: 'FIN', prop: 'Allele_frequency_exome_FIN_GnomAD'},
                {label: 'NFE', prop: 'Allele_frequency_exome_NFE_GnomAD'},
                {label: 'SAS', prop: 'Allele_frequency_exome_SAS_GnomAD'},
                {label: 'OTH', prop: 'Allele_frequency_exome_OTH_GnomAD'},
            ];
            title = 'gnomAD Exomes';
            pointFormat =  "{point.y}<br /><em>({point.count} of {point.number})</em>";
        }
        else if (prop === 'Allele_Frequency_Charts_Genome_GnomADv3') {
            if (!variant.Variant_in_GnomADv3) {
                return false;
            }
            frequencyProps = [
                {label: 'AFR', prop: 'Allele_frequency_genome_AFR_GnomADv3'},
                {label: 'AMR', prop: 'Allele_frequency_genome_AMR_GnomADv3'},
                {label: 'AMI', prop: 'Allele_frequency_genome_AMI_GnomADv3'},
                {label: 'ASJ', prop: 'Allele_frequency_genome_ASJ_GnomADv3'},
                {label: 'EAS', prop: 'Allele_frequency_genome_EAS_GnomADv3'},
                {label: 'FIN', prop: 'Allele_frequency_genome_FIN_GnomADv3'},
                {label: 'MID', prop: 'Allele_frequency_genome_MID_GnomADv3'},
                {label: 'NFE', prop: 'Allele_frequency_genome_NFE_GnomADv3'},
                {label: 'SAS', prop: 'Allele_frequency_genome_SAS_GnomADv3'},
                {label: 'OTH', prop: 'Allele_frequency_genome_OTH_GnomADv3'},
            ];
            title = 'gnomAD Genomes';
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
            <tr>
                <td>
                    <div className="alleleFrequencyChartOuterContainer">
                        <div className="alleleFrequencyChartContainer">
                            <div className='alleleFrequencyChart'>
                                <BarChart ref={`${prop}_alleleFreq2`} container={`${prop}_alleleFreq1`} options={fullscaleChartOptions}/>
                            </div>
                            <div className='alleleFrequencyChart' onClick={this.toggleScale.bind(this, `${prop}_alleleFreq2`)}>
                                <BarChart ref={`${prop}_alleleFreq2`} container={`${prop}_alleleFreq2`} options={scaledChartOptions} />
                                <div style={{textAlign: 'center', color: 'grey', fontSize: '12px'}}>(click chart to change scale)</div>
                            </div>
                        </div>
                    </div>
                </td>
            </tr>
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
